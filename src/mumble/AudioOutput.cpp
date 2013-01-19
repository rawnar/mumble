/* Copyright (C) 2005-2011, Thorvald Natvig <thorvald@natvig.com>
   Copyright (C) 2009-2011, Stefan Hacker <dd0t@users.sourceforge.net>

   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   - Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.
   - Neither the name of the Mumble Developers nor the names of its
     contributors may be used to endorse or promote products derived from this
     software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "mumble_pch.hpp"

#include "AudioOutput.h"

#include "AudioOutputSample.h"
#include "AudioOutputSpeech.h"
#include "User.h"
#include "Global.h"
#include "Message.h"
#include "Plugins.h"
#include "PacketDataStream.h"
#include "ServerHandler.h"
#include "VoiceRecorder.h"

// Remember that we cannot use static member classes that are not pointers, as the constructor
// for AudioOutputRegistrar() might be called before they are initialized, as the constructor
// is called from global initialization.
// Hence, we allocate upon first call.

QMap<QString, AudioOutputRegistrar *> *AudioOutputRegistrar::qmNew;
QString AudioOutputRegistrar::current = QString();

AudioOutputRegistrar::AudioOutputRegistrar(const QString &n, int p) : name(n), priority(p) {
	if (! qmNew)
		qmNew = new QMap<QString, AudioOutputRegistrar *>();
	qmNew->insert(name,this);
}

AudioOutputRegistrar::~AudioOutputRegistrar() {
	qmNew->remove(name);
}

AudioOutputPtr AudioOutputRegistrar::newFromChoice(QString choice) {
	if (! qmNew)
		return AudioOutputPtr();

	if (!choice.isEmpty() && qmNew->contains(choice)) {
		g.s.qsAudioOutput = choice;
		current = choice;
		return AudioOutputPtr(qmNew->value(choice)->create());
	}
	choice = g.s.qsAudioOutput;
	if (qmNew->contains(choice)) {
		current = choice;
		return AudioOutputPtr(qmNew->value(choice)->create());
	}

	AudioOutputRegistrar *r = NULL;
	foreach(AudioOutputRegistrar *aor, *qmNew)
		if (!r || (aor->priority > r->priority))
			r = aor;
	if (r) {
		current = r->name;
		return AudioOutputPtr(r->create());
	}
	return AudioOutputPtr();
}

bool AudioOutputRegistrar::canMuteOthers() const {
	return false;
}

bool AudioOutputRegistrar::usesOutputDelay() const {
	return true;
}

bool AudioOutputRegistrar::canExclusive() const {
	return false;
}

AudioOutput::AudioOutput() {
	iFrameSize = SAMPLE_RATE / 100;
	bRunning = true;

	iChannels = 0;
	fSpeakers = NULL;
	bSpeakerPositional = NULL;
	bSpeakerDirectional = NULL;
	iSpeakerMap = NULL;

	iMixerFreq = 0;
	eSampleFormat = SampleFloat;
	iSampleSize = 0;
}

AudioOutput::~AudioOutput() {
	bRunning = false;
	wait();
	wipe();

	delete [] fSpeakers;
	delete [] bSpeakerPositional;
	delete [] bSpeakerDirectional;
	delete [] iSpeakerMap;
}

void AudioOutput::wipe() {
	foreach(AudioOutputUser *aop, qmOutputs)
		removeBuffer(aop);
}

const float *AudioOutput::getSpeakerPos(unsigned int &speakers) {
	if ((iChannels > 0) && fSpeakers) {
		speakers = iChannels;
		return fSpeakers;
	}
	return NULL;
}

AudioOutputSample *AudioOutput::playSample(const QString &filename, bool loop) {
	SoundFile *handle = AudioOutputSample::loadSndfile(filename);
	if (handle == NULL)
		return NULL;

	while ((iMixerFreq == 0) && isAlive()) {
		QThread::yieldCurrentThread();
	}

	if (! iMixerFreq)
		return NULL;

	QWriteLocker locker(&qrwlOutputs);
	AudioOutputSample *aos = new AudioOutputSample(filename, handle, loop, iMixerFreq);
	qmOutputs.insert(NULL, aos);

	return aos;

}

void AudioOutput::addFrameToBuffer(ClientUser *user, const QByteArray &qbaPacket, unsigned int iSeq, MessageHandler::UDPMessageType type) {
	if (iChannels == 0)
		return;
	qrwlOutputs.lockForRead();
	AudioOutputSpeech *aop = qobject_cast<AudioOutputSpeech *>(qmOutputs.value(user));

	if (! aop || (aop->umtType != type)) {
		qrwlOutputs.unlock();

		if (aop)
			removeBuffer(aop);

		while ((iMixerFreq == 0) && isAlive()) {
			QThread::yieldCurrentThread();
		}

		if (! iMixerFreq)
			return;

		qrwlOutputs.lockForWrite();
		aop = new AudioOutputSpeech(user, iMixerFreq, type);
		qmOutputs.replace(user, aop);
	}

	aop->addFrameToBuffer(qbaPacket, iSeq);

	qrwlOutputs.unlock();
}

void AudioOutput::removeBuffer(const ClientUser *user) {
	removeBuffer(qmOutputs.value(user));
}

void AudioOutput::removeBuffer(AudioOutputUser *aop) {
	QWriteLocker locker(&qrwlOutputs);
	for (QMultiHash<const ClientUser *, AudioOutputUser *>::iterator i=qmOutputs.begin();
			i!=qmOutputs.end();++i) {
		if (i.value() == aop) {
			qmOutputs.erase(i);
			delete aop;
			break;
		}
	}
}

/*! \brief A function that takes the chanmasks and positions the loudspeakers around the user.
 *
 * This function will position the loudspeakers around the user. It will also distinguish between
 * loudspeakers that will not be used for the positional audio and loudspeakers that have no
 * directional value like the LFE. Speakers that are not in the xz-plane are not yet used.
 *
 * \param chanmasks a constant unsigned int holding the channel masks
 * \param forceheadphone a boolean that will force the use of headphones
 */
void AudioOutput::setSpeakerPositions(const unsigned int *chanmasks, bool forceheadphone) {
	delete [] fSpeakers;
	delete [] bSpeakerPositional;
	delete [] bSpeakerDirectional;
	delete [] iSpeakerMap;

	fSpeakers = new float[iChannels * 3];
	bSpeakerPositional = new bool[iChannels];
	bSpeakerDirectional = new bool[iChannels];

	memset(fSpeakers, 0, sizeof(float) * iChannels * 3);
	memset(bSpeakerPositional, 0, sizeof(bool) * iChannels);
	memset(bSpeakerDirectional, 0, sizeof(bool) * iChannels);

	if (g.s.bPositionalHeadphone || forceheadphone) {
		bHeadphones = true;
	} else {
		bHeadphones = false;
	}

	// This QMap will be used to sort the loudspeakers counter-clockwise. The loudspeakers
	// are numbered counter-clockwise starting with the side right loudspeaker.
	QMap<int, unsigned int> map;
	for (unsigned int i=0;i<iChannels;++i) {
		float *s = &fSpeakers[3*i];
		bSpeakerPositional[i] = true;
		bSpeakerDirectional[i] = true;

		switch (chanmasks[i]) {
			case SPEAKER_FRONT_LEFT:
				s[0] = -0.5f;
				s[2] = 1.0f;
				map[6] = i;
				break;
			case SPEAKER_FRONT_RIGHT:
				s[0] = 0.5f;
				s[2] = 1.0f;
				map[2] = i;
				break;
			case SPEAKER_FRONT_CENTER:
				s[2] = 1.0f;
				map[4] = i;
				break;
			case SPEAKER_LOW_FREQUENCY:
				bSpeakerDirectional[i] = false;
				break;
			case SPEAKER_BACK_LEFT:
				s[0] = -0.5f;
				s[2] = -1.0f;
				map[8] = i;
				break;
			case SPEAKER_BACK_RIGHT:
				s[0] = 0.5f;
				s[2] = -1.0f;
				map[10] = i;
				break;
			case SPEAKER_FRONT_LEFT_OF_CENTER:
				s[0] = -0.25;
				s[2] = 1.0f;
				map[5] = i;
				break;
			case SPEAKER_FRONT_RIGHT_OF_CENTER:
				s[0] = 0.25;
				s[2] = 1.0f;
				map[3] = i;
				break;
			case SPEAKER_BACK_CENTER:
				s[2] = -1.0f;
				map[9] = i;
				break;
			case SPEAKER_SIDE_LEFT:
				s[0] = -1.0f;
				map[7] = i;
				break;
			case SPEAKER_SIDE_RIGHT:
				s[0] = 1.0f;
				map[1] = i;
				break;
			case SPEAKER_TOP_CENTER:
			case SPEAKER_TOP_FRONT_LEFT:
			case SPEAKER_TOP_FRONT_CENTER:
			case SPEAKER_TOP_FRONT_RIGHT:
			case SPEAKER_TOP_BACK_LEFT:
			case SPEAKER_TOP_BACK_CENTER:
			case SPEAKER_TOP_BACK_RIGHT:
				bSpeakerPositional[i] = false;
				qWarning("AudioOutput: Out of plane speakers are not yet supported");
				break;
			default:
				bSpeakerPositional[i] = false;
				qWarning("AudioOutput: Unknown speaker %d: %08x", i, chanmasks[i]);
				break;
		}
		if (bHeadphones) {
			s[1] = 0.0f;
			s[2] = 0.0f;
			if (s[0] == 0.0f && bSpeakerDirectional[i])
				bSpeakerPositional[i] = false;
		}
	}
	// The amount of positional speakers
	iChannelsPositional = map.size();
	iSpeakerMap = new unsigned int[iChannelsPositional + 2];
	// The keys are sorted in a QMap when going through the iterator
	unsigned int m = 0;
	for (QMap<int, unsigned int>::const_iterator i=map.constBegin();i!=map.constEnd();++i) {
		++m;
		iSpeakerMap[m] = i.value();
		//qWarning("m = %d,key = %d, value = %d",m,i.key(),i.value());
	}
	// Closing the mapping, so it will be easier to use
	iSpeakerMap[0] = iSpeakerMap[m];
	iSpeakerMap[m+1] = iSpeakerMap[1];

	// Normalizing the loudspeaker distance
	for (unsigned int i=0;i<iChannels;++i) {
		float d = sqrtf(fSpeakers[3*i+0]*fSpeakers[3*i+0] + fSpeakers[3*i+1]*fSpeakers[3*i+1] + fSpeakers[3*i+2]*fSpeakers[3*i+2]);
		if (d > 0.0f) {
			fSpeakers[3*i+0] /= d;
			fSpeakers[3*i+1] /= d;
			fSpeakers[3*i+2] /= d;
		}
	}
}

void AudioOutput::initializeMixer(const unsigned int *chanmasks, bool forceheadphone) {
	if (g.s.bPositionalAudio && (iChannels > 1))
		setSpeakerPositions(chanmasks, forceheadphone);

	iSampleSize = static_cast<int>(iChannels * ((eSampleFormat == SampleFloat) ? sizeof(float) : sizeof(short)));
	qWarning("AudioOutput: Initialized %d channel %d hz mixer", iChannels, iMixerFreq);
}

/*! \brief Function that will postion the listener in a 3-D world.
 *
 * The listener is positioned at its location and the local coordinate frame of his head is build.
 *
 * \param position this is a vector that gives the listeners location in the global frame
 * \param front    the front vector of the head given in the global frame
 * \param top      the top vector of the head given in the global frame
 *
 * \return boolean depending on success
 */
bool AudioOutput::positionListener(float *position, float *front, float *top) {

	fPosition[0] = position[0];
	fPosition[1] = position[1];
	fPosition[2] = position[2];

	// Front vector is dominant; if it's zero we presume all is zero

	float len = sqrtf(front[0]*front[0]+front[1]*front[1]+front[2]*front[2]);

	if (len > 0.0f) {
		fFront[0] = front[0] / len;
		fFront[1] = front[1] / len;
		fFront[2] = front[2] / len;

		len = sqrtf(top[0]*top[0]+top[1]*top[1]+top[2]*top[2]);

		if (len > 0.0f) {
			fTop[0] = top[0] / len;
			fTop[1] = top[1] / len;
			fTop[2] = top[2] / len;
		} else {
			fTop[0] = 0.0f;
			fTop[1] = 1.0f;
			fTop[2] = 0.0f;
		}

		if (fabs(fFront[0] * fTop[0] + fFront[1] * fTop[1] + fFront[2] * fTop[2]) > 0.01f) {
			// Not perpendicular. Assume Y up and rotate 90 degrees
			float azimuth = 0.0f;
			if ((fabs(front[0]) >= 0.0f) || (fabs(front[2]) >= 0.0f))
				azimuth = atan2f(front[2], front[0]);
			float inclination = acosf(front[1]) - float(M_PI) / 2.0f;

			fTop[0] = sinf(inclination) * cosf(azimuth);
			fTop[1] = cosf(inclination);
			fTop[2] = sinf(inclination) * sinf(azimuth);
		}
		// Calculate right vector as front X top
		fRight[0] = fTop[1] * fFront[2] - fTop[2] * fFront[1];
		fRight[1] = fTop[2] * fFront[0] - fTop[0] * fFront[2];
		fRight[2] = fTop[0] * fFront[1] - fTop[1] * fFront[0];

	} else {
		fRight[0] = 1.0f;
		fRight[1] = 0.0f;
		fRight[2] = 0.0f;

		fTop[0] = 0.0f;
		fTop[1] = 1.0f;
		fTop[2] = 0.0f;

		fFront[0] = 0.0f;
		fFront[1] = 0.0f;
		fFront[2] = 1.0f;
	}

	/*
	qWarning("Front: %f %f %f", fFront[0], fFront[1], fFront[2]);
	qWarning("Top: %f %f %f", fTop[0], fTop[1], fTop[2]);
	qWarning("Right: %f %f %f", fRight[0], fRight[1], fRight[2]);
	*/

	return true;
}

/*! \brief sets the location of the audio source.
 *
 * This method will set the location of the audio source and calculates the direction and distance
 * relative to the listener. It will also calculate the distance attenuation and the bloom
 * amplification needed later on.
 *
 * \param position a vector that gives the location of the audio source in the global frame
 *
 * \return boolean will be true when it is a viable source, false if not.
 */
boost::optional<AudioSourceData> AudioOutput::locateSource(const float * const position) {
	// calculate position relative to own position
	float relSrcPosition[3] = { position[0] - fPosition[0], position[1] - fPosition[1], position[2] - fPosition[2] };
	const float distance = sqrtf(relSrcPosition[0] * relSrcPosition[0] + relSrcPosition[1] * relSrcPosition[1] + relSrcPosition[2] * relSrcPosition[2]);
	if (distance > 0.01f) {
		relSrcPosition[0] /= distance;
		relSrcPosition[1] /= distance;
		relSrcPosition[2] /= distance;
	} else {
		// Source is very close to the listener and thus the source is not viable
		boost::optional<AudioSourceData> retData;
		return retData;
	}

	boost::optional<AudioSourceData> retData;
	retData = AudioSourceData();

	// rotate the direction vector into the local frame of the listener
	retData->fDirection[0] = relSrcPosition[0] * fRight[0] + relSrcPosition[1] * fRight[1] + relSrcPosition[2] * fRight[2];
	retData->fDirection[1] = relSrcPosition[0] * fTop[0]   + relSrcPosition[1] * fTop[1]   + relSrcPosition[2] * fTop[2];
	retData->fDirection[2] = relSrcPosition[0] * fFront[0] + relSrcPosition[1] * fFront[1] + relSrcPosition[2] * fFront[2];

	// distance attenuation
	if (g.s.fAudioMaxDistVolume > 0.99f || distance < g.s.fAudioMinDistance) {
		retData->fAttenuation = 1.0f;
	} else if (distance >= g.s.fAudioMaxDistance) {
		retData->fAttenuation = g.s.fAudioMaxDistVolume;
	} else {
		float mvol = g.s.fAudioMaxDistVolume;
		if (mvol < 0.01f)
			mvol = 0.01f;

		const float drel = (distance - g.s.fAudioMinDistance) / (g.s.fAudioMaxDistance - g.s.fAudioMinDistance);
		retData->fAttenuation = powf(10.0f, log10f(mvol) * drel);
	}

	// Here's the theory.
	// We support sound "bloom"ing. That is, if sound comes directly from the left, if it is sufficiently
	// close, we'll hear it full intensity from the left side, and "bloom" intensity from the right side.
	if (distance < g.s.fAudioMinDistance) {
		retData->fBloom = g.s.fAudioBloom * (1.0f - distance/g.s.fAudioMinDistance);
	} else {
		retData->fBloom = 0.0f;
	}

	return retData;
}

/*! \brief Calculates the gain of a given loudspeaker depending on the direction of the sound
 *
 * It calculates the gain of the given loudspeaker. The audio will only be distributed to the
 * loudspeakers that will clamp the direction vector.
 *
 * \param chanIndex the number of the loudspeaker
 * \return the gain of the loudspeaker
 * \note Elevation is not taken into account, as it is non distinguishable from distance.
 */
float AudioOutput::calcGain(const unsigned int & chanIndex, const AudioSourceData & sourceData) {
	// If speakers are not positional remove them
	if (!bSpeakerPositional[chanIndex]) {
		return 0.0f;
	}

	// If speakers are non directional only use distance attenuation
	if (!bSpeakerDirectional[chanIndex]) {
		return qMin(1.0f, sourceData.fAttenuation + sourceData.fBloom);
	}

	// For headphones fall-back to the old method
	if (bHeadphones) {
		float dotProduct = sourceData.fDirection[0] * fSpeakers[chanIndex*3+0]
						   + sourceData.fDirection[1] * fSpeakers[chanIndex*3+1]
						   + sourceData.fDirection[2] * fSpeakers[chanIndex*3+2];
		float dotFactor = (dotProduct + 1.0f) / 2.0f;
		return qMin(1.0f, sourceData.fAttenuation * dotFactor + sourceData.fBloom);
	}

	// Rotate the direction vector to the xz-plane
	const float len = sqrtf(
			sourceData.fDirection[0] * sourceData.fDirection[0]
			+ sourceData.fDirection[2] * sourceData.fDirection[2]);

	// This happens when fDirection is pointing along the y-axis
	if (len <= 0.01f) {
		return qMin(1.0f, sourceData.fAttenuation / float(iChannelsPositional) + sourceData.fBloom);
	}

	const float xp = sourceData.fDirection[0] / len;
	const float zp = sourceData.fDirection[2] / len;

	// Find the loudspeakers right and left of the loudspeaker in question
	const float * const sc = &fSpeakers[3*chanIndex];
	float * sr = NULL;
	float * sl = NULL;
	for (unsigned int i=1;i<iChannelsPositional+1;++i) {
		if (iSpeakerMap[i] == chanIndex) {
			sr = &fSpeakers[3*iSpeakerMap[i-1]];
			sl = &fSpeakers[3*iSpeakerMap[i+1]];
			break;
		}
	}
	if (!sr && !sl) // loudspeaker was not found, return zero
		return 0.0f;

	// Find out if the source is located between sc-sr or sc-sl
	float dot12;
	float dot1;
	float dot2;
	if ((sr[0] - xp) * (sc[2] - zp) <= (sr[2] - zp) * (sc[0] - xp)) { // point p is on the rightside of line sr to sc
		dot12 = sc[0] * sr[0] + sc[2] * sr[2];
		dot1 = sc[0] * xp + sc[2] * zp;
		dot2 = sr[0] * xp + sr[2] * zp;
	} else if ((sc[0] - xp) * (sl[2] - zp) <= (sc[2] - zp) * (sl[0] - xp)) { // point p is on the rightside of line sc to sl
		dot12 = sc[0] * sl[0] + sc[2] * sl[2];
		dot1 = sc[0] * xp + sc[2] * zp;
		dot2 = sl[0] * xp + sl[2] * zp;
	} else {
		return qMin(1.0f, sourceData.fBloom);
	}

	// Calculate the gain for the loudspeaker in question
	if (dot1 + dot2 >= 1.0f + dot12) { // normal case angle between loudspeakers is smaller then 180 degrees
		float dotFactor = acosf(dot2) / (acosf(dot1) + acosf(dot2));
		return qMin(1.0f, sourceData.fAttenuation * dotFactor + sourceData.fBloom);
	} else if (dot12 + dot1 + dot2 <= -1.0f) { // source is between loudspeakers but on the opposite side of the listener
		float dotFactor = acosf(-dot1) / (acosf(-dot1) + acosf(-dot2));
		return qMin(1.0f, sourceData.fAttenuation * dotFactor + sourceData.fBloom);
	} else if ( dot1 > dot2 ) { // special case were all power goes to the loudspeaker in question
		return qMin(1.0f, sourceData.fAttenuation + sourceData.fBloom);
	} else { // when all else fails only bloom is needed
		return qMin(1.0f, sourceData.fBloom);
	}
}

bool AudioOutput::mix(void *outbuff, unsigned int nsamp) {
	QList<AudioOutputUser *> qlMix;
	QList<AudioOutputUser *> qlDel;

	if (g.s.fVolume < 0.01f)
		return false;

	const float adjustFactor = std::pow(10, -18.f / 20);
	const float defaultVolume = g.s.fVolume;
	const unsigned int nchan = iChannels;
	ServerHandlerPtr sh = g.sh;
	VoiceRecorderPtr recorder;
	if (sh)
		recorder = g.sh->recorder;

	qrwlOutputs.lockForRead();
	bool needAdjustment = false;
	for (QMultiHash<const ClientUser *, AudioOutputUser *>::const_iterator it=qmOutputs.constBegin();
		       it!=qmOutputs.constEnd();++it) {
		AudioOutputUser *aop = it.value();
		if (! aop->needSamples(nsamp)) {
			qlDel.append(aop);
		} else {
			qlMix.append(aop);
			// Set a flag if there is a priority speaker
			if (it.key() && it.key()->bPrioritySpeaker)
				needAdjustment = true;
		}
	}

	if (! qlMix.isEmpty()) {
		STACKVAR(float, fOutput, iChannels * nsamp);
		float *output = (eSampleFormat == SampleFloat) ? reinterpret_cast<float *>(outbuff) : fOutput;
		bool validListener = false;

		memset(output, 0, sizeof(float) * nsamp * iChannels);

		boost::shared_array<float> recbuff;
		if (recorder) {
			recbuff = boost::shared_array<float>(new float[nsamp]);
			memset(recbuff.get(), 0, sizeof(float) * nsamp);
		}

		if (g.s.bPositionalAudio && (iChannels > 1) && g.p->fetch() && (g.bPosTest || fabs(g.p->fCameraPosition[0]) > 0.0f || fabs(g.p->fCameraPosition[1]) > 0.0f || fabs(g.p->fCameraPosition[2]) > 0.0f))
			validListener = positionListener(g.p->fCameraPosition, g.p->fCameraFront, g.p->fCameraTop);

		foreach(AudioOutputUser *aop, qlMix) {
			const float * RESTRICT pfBuffer = aop->pfBuffer;
			float prioSpeakerVolAdjFactor = 1.0; // holds volume adjustment in case of priority speaker
			boost::optional<AudioSourceData> srcData = NULL;

			// If we have at least one priority speaker
			if (needAdjustment) {
				AudioOutputSpeech *aos = qobject_cast<AudioOutputSpeech *>(aop);
				// Exclude whispering people
				if (aos && (aos->p->tsState == Settings::Talking || aos->p->tsState == Settings::Shouting)) {
					// Adjust all non-priority speakers
					if (!aos->p->bPrioritySpeaker)
						prioSpeakerVolAdjFactor = adjustFactor;
				}
			}

			if (recorder) {
				AudioOutputSpeech *aos = qobject_cast<AudioOutputSpeech *>(aop);

				if (aos) {
					for (unsigned int i=0;i<nsamp;++i) {
						recbuff[i] += pfBuffer[i] * prioSpeakerVolAdjFactor;
					}

					if (!recorder->getMixDown()) {
						if (aos) {
							recorder->addBuffer(aos->p, recbuff, nsamp);
						} else {
							// this should be unreachable
							Q_ASSERT(false);
						}
						recbuff = boost::shared_array<float>(new float[nsamp]);
						memset(recbuff.get(), 0, sizeof(float) * nsamp);
					}

					// Don't add the local audio to the real output
					if (qobject_cast<RecordUser *>(aos->p)) {
						continue;
					}
				}
			}

			if (validListener && (fabs(aop->fPos[0]) > 0.0f || fabs(aop->fPos[1]) > 0.0f || fabs(aop->fPos[2]) > 0.0f))
				srcData = locateSource(aop->fPos);

			if (srcData) {
				// aop->fPos is a viable source for pos audio
				// if pfVolume is uninitialized, allocate for number of channels and initialize values with -1
				if (! aop->pfVolume) {
					aop->pfVolume = new float[nchan];
					for (unsigned int s=0;s<nchan;++s)
						aop->pfVolume[s] = -1.0;
				}
				for (unsigned int chanIndex=0; chanIndex < nchan; ++chanIndex) {
					const float volumeAdjFactor = defaultVolume * calcGain(chanIndex, *srcData) * prioSpeakerVolAdjFactor;
					float * RESTRICT o = output + chanIndex;
					// if pfVolume is unset, use volumeAdjFactor as old, otherwise the pfVolume value
					const float old = (aop->pfVolume[chanIndex] >= 0.0) ? aop->pfVolume[chanIndex] : volumeAdjFactor;
					const float inc = (volumeAdjFactor - old) / static_cast<float>(nsamp);
					aop->pfVolume[chanIndex] = volumeAdjFactor;
					// pass data to output buffer
					if ((old >= 0.00000001f) || (volumeAdjFactor >= 0.00000001f))
						for (unsigned int i=0;i<nsamp;++i)
							o[i*nchan] += pfBuffer[i] * (old + inc * static_cast<float>(i));
				}
			} else {
				// Just plain old 1D audio
				// For each channel adjust buffered data (from output user) with volume factor and set in output buffer
				for (unsigned int s=0;s<nchan;++s) {
					const float volumeAdjFactor = defaultVolume * prioSpeakerVolAdjFactor;
					float * RESTRICT o = output + s;
					for (unsigned int i=0;i<nsamp;++i)
						o[i*nchan] += pfBuffer[i] * volumeAdjFactor;
				}
			}
		}

		if (recorder && recorder->getMixDown()) {
			recorder->addBuffer(NULL, recbuff, nsamp);
		}

		// Clip
		if (eSampleFormat == SampleFloat)
			for (unsigned int i=0;i<nsamp*iChannels;++i)
				output[i] = qBound(-1.0f, output[i], 1.0f);
		else
			for (unsigned int i=0;i<nsamp*iChannels;++i)
				reinterpret_cast<short *>(outbuff)[i] = static_cast<short>(qBound(-32768.f, (output[i] * 32768.f), 32767.f));
	}

	qrwlOutputs.unlock();

	foreach(AudioOutputUser *aop, qlDel)
		removeBuffer(aop);

	return (! qlMix.isEmpty());
}

bool AudioOutput::isAlive() const {
	return isRunning();
}

unsigned int AudioOutput::getMixerFreq() const {
	return iMixerFreq;
}
