//
//  Reverb.hpp
//  JuceDemoPlugin
//
//  Created by yuasa on 2017/01/23.
//
//

#ifndef Reverb_hpp
#define Reverb_hpp

#include <stdio.h>

#include "../JuceLibraryCode/JuceHeader.h"

void SetTimeParam(float seconds);
void SetErTime(float milliseconds);
void SetAirAbsorption(float freq);
void SetStereoSpread(float amount);
void SetErDamping(float frequency);
void SetPreDelayTime(float milliseconds);
void SetFeedbackGain(float feedback_gain);

void ApplyReverb(AudioBuffer<float> const & input,
                 AudioBuffer<float> & output,
                 int length,
                 float dry_wet,
                 float er_gain);

#endif /* Reverb_hpp */
