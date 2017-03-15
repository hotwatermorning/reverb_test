/*
  ==============================================================================

    This file was auto-generated by the Jucer!

    It contains the basic startup code for a Juce application.

  ==============================================================================
*/

#ifndef __PLUGINEDITOR_H_4ACCBAA__
#define __PLUGINEDITOR_H_4ACCBAA__

#include "../JuceLibraryCode/JuceHeader.h"
#include "PluginProcessor.h"


//==============================================================================
/** This is the editor component that our filter will display.
*/
class ReverbTestAudioProcessorEditor  : public AudioProcessorEditor,
                                            private Timer
{
public:
    ReverbTestAudioProcessorEditor (ReverbTestAudioProcessor&);
    ~ReverbTestAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;
    void timerCallback() override;

private:
    class ParameterSlider;

    Label timecodeDisplayLabel, buildTimeLabel, wetLabel, timeLabel, lpfLabel,
    erTimeLabel, erGainLabel, erDampingLabel,
    stereoSpreadLabel, preDelayTimeLabel,
    feedbackGainLabel;
    
    ScopedPointer<ParameterSlider> wetSlider, timeSlider, lpfSlider,
    erTimeSlider, erGainSlider, erDampingSlider,
    stereoSpreadSlider, preDelayTimeSlider,
    feedbackGainSlider;
    
    ScopedPointer<ResizableCornerComponent> resizer;
    ComponentBoundsConstrainer resizeLimits;

    AudioPlayHead::CurrentPositionInfo lastDisplayedPosition;

    //==============================================================================
    ReverbTestAudioProcessor& getProcessor() const
    {
        return static_cast<ReverbTestAudioProcessor&> (processor);
    }

    void updateTimecodeDisplay (AudioPlayHead::CurrentPositionInfo);
};


#endif  // __PLUGINEDITOR_H_4ACCBAA__
