/*
  ==============================================================================

    This file was auto-generated by the Jucer!

    It contains the basic startup code for a Juce application.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
// This is a handy slider subclass that controls an AudioProcessorParameter
// (may move this class into the library itself at some point in the future..)
class ReverbTestAudioProcessorEditor::ParameterSlider
:   public Slider
,   private Timer
{
public:
    ParameterSlider (AudioProcessorParameter& p)
        : Slider (p.getName (256)), param (p)
    {
        setRange (0.0, 1.0, 0.0);
        startTimerHz (30);
        updateSliderPos();
    }

    void valueChanged() override
    {
        param.setValue ((float) Slider::getValue());
    }

    void timerCallback() override       { updateSliderPos(); }

    void startedDragging() override     { param.beginChangeGesture(); }
    void stoppedDragging() override     { param.endChangeGesture();   }

    double getValueFromText (const String& text) override   { return param.getValueForText (text); }
    String getTextFromValue (double value) override         { return param.getText ((float) value, 1024); }

    void updateSliderPos()
    {
        const float newValue = param.getValue();

        if (newValue != (float) Slider::getValue())
            Slider::setValue (newValue);
    }

    AudioProcessorParameter& param;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ParameterSlider)
};

//==============================================================================
ReverbTestAudioProcessorEditor::ReverbTestAudioProcessorEditor (ReverbTestAudioProcessor& owner)
    : AudioProcessorEditor (owner),
      timecodeDisplayLabel (String::empty),
      wetLabel (String::empty, "Dry/Wet:"),
      timeLabel (String::empty, "Time:"),
      lpfLabel (String::empty, "Air Absorption"),
      erTimeLabel(String::empty, "ER Time"),
      erGainLabel(String::empty, "ER Gain"),
      erDampingLabel(String::empty, "Er Damping"),
      stereoSpreadLabel(String::empty, "Stereo Spread"),
      preDelayTimeLabel(String::empty, "Pre Delay Time(ms)")
{
    // add some sliders..
    addAndMakeVisible (wetSlider = new ParameterSlider (*owner.wetParam));
    wetSlider->setSliderStyle (Slider::Rotary);

    addAndMakeVisible (timeSlider = new ParameterSlider (*owner.timeParam));
    timeSlider->setSliderStyle (Slider::Rotary);
    
    addAndMakeVisible (lpfSlider = new ParameterSlider (*owner.lpfParam));
    lpfSlider->setSliderStyle (Slider::Rotary);

    addAndMakeVisible (erTimeSlider = new ParameterSlider (*owner.erTimeParam));
    erTimeSlider->setSliderStyle (Slider::Rotary);
    
    addAndMakeVisible (erGainSlider = new ParameterSlider (*owner.erGainParam));
    erGainSlider->setSliderStyle (Slider::Rotary);
    
    addAndMakeVisible (erDampingSlider = new ParameterSlider (*owner.erDampingParam));
    erDampingSlider->setSliderStyle (Slider::Rotary);
    
    addAndMakeVisible (stereoSpreadSlider = new ParameterSlider (*owner.stereoSpreadParam));
    stereoSpreadSlider->setSliderStyle (Slider::Rotary);
    
    addAndMakeVisible (preDelayTimeSlider = new ParameterSlider (*owner.preDelayParam));
    preDelayTimeSlider->setSliderStyle (Slider::Rotary);

    // add some labels for the sliders..
    wetLabel.attachToComponent (wetSlider, false);
    wetLabel.setFont (Font (11.0f));

    timeLabel.attachToComponent (timeSlider, false);
    timeLabel.setFont (Font (11.0f));
    
    lpfLabel.attachToComponent (lpfSlider, false);
    lpfLabel.setFont (Font (11.0f));

    erTimeLabel.attachToComponent (erTimeSlider, false);
    erTimeLabel.setFont (Font (11.0f));
    
    erGainLabel.attachToComponent (erGainSlider, false);
    erGainLabel.setFont (Font (11.0f));
    
    erDampingLabel.attachToComponent (erDampingSlider, false);
    erDampingLabel.setFont (Font (11.0f));
    
    stereoSpreadLabel.attachToComponent (stereoSpreadSlider, false);
    stereoSpreadLabel.setFont (Font (11.0f));
    
    preDelayTimeLabel.attachToComponent (preDelayTimeSlider, false);
    preDelayTimeLabel.setFont (Font (11.0f));

    // add a label that will display the current timecode and status..
    addAndMakeVisible (timecodeDisplayLabel);
    timecodeDisplayLabel.setColour (Label::textColourId, Colours::blue);
    timecodeDisplayLabel.setFont (Font (Font::getDefaultMonospacedFontName(), 15.0f, Font::plain));

    // add the triangular resizer component for the bottom-right of the UI
    addAndMakeVisible (resizer = new ResizableCornerComponent (this, &resizeLimits));
    resizeLimits.setSizeLimits (150, 150, 800, 300);

    // set our component's initial size to be the last one that was stored in the filter's settings
    setSize (owner.lastUIWidth,
             owner.lastUIHeight);

    // start a timer which will keep our timecode display updated
    startTimerHz (30);
}

ReverbTestAudioProcessorEditor::~ReverbTestAudioProcessorEditor()
{
}

//==============================================================================
void ReverbTestAudioProcessorEditor::paint (Graphics& g)
{
    g.setGradientFill (ColourGradient (Colours::white, 0, 0,
                                       Colours::lightgrey, 0, (float) getHeight(), false));
    g.fillAll();
}

void ReverbTestAudioProcessorEditor::resized()
{
    // This lays out our child components...

    Rectangle<int> r (getLocalBounds().reduced (8));

    timecodeDisplayLabel.setBounds (r.removeFromTop (26));

    r.removeFromTop (30);
    Rectangle<int> sliderArea (r.removeFromTop (50));
    wetSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 8)));
    timeSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 7)));
    lpfSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 6)));
    erTimeSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 5)));
    erGainSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 4)));
    erDampingSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 3)));
    stereoSpreadSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 2)));
    preDelayTimeSlider->setBounds (sliderArea.removeFromLeft (jmin (180, sliderArea.getWidth() / 1)));

    resizer->setBounds (getWidth() - 16, getHeight() - 16, 16, 16);

    getProcessor().lastUIWidth = getWidth();
    getProcessor().lastUIHeight = getHeight();
}

//==============================================================================
void ReverbTestAudioProcessorEditor::timerCallback()
{
    updateTimecodeDisplay (getProcessor().lastPosInfo);
}

//==============================================================================
// quick-and-dirty function to format a timecode string
static String timeToTimecodeString (double seconds)
{
    const int millisecs = roundToInt (std::abs (seconds * 1000.0));

    return String::formatted ("%s%02d:%02d:%02d.%03d",
                              seconds < 0 ? "-" : "",
                              millisecs / 360000,
                              (millisecs / 60000) % 60,
                              (millisecs / 1000) % 60,
                              millisecs % 1000);
}

// quick-and-dirty function to format a bars/beats string
static String quarterNotePositionToBarsBeatsString (double quarterNotes, int numerator, int denominator)
{
    if (numerator == 0 || denominator == 0)
        return "1|1|000";

    const int quarterNotesPerBar = (numerator * 4 / denominator);
    const double beats  = (fmod (quarterNotes, quarterNotesPerBar) / quarterNotesPerBar) * numerator;

    const int bar    = ((int) quarterNotes) / quarterNotesPerBar + 1;
    const int beat   = ((int) beats) + 1;
    const int ticks  = ((int) (fmod (beats, 1.0) * 960.0 + 0.5));

    return String::formatted ("%d|%d|%03d", bar, beat, ticks);
}

// Updates the text in our position label.
void ReverbTestAudioProcessorEditor::updateTimecodeDisplay (AudioPlayHead::CurrentPositionInfo pos)
{
    if (lastDisplayedPosition != pos)
    {
        lastDisplayedPosition = pos;

        MemoryOutputStream displayText;

        displayText << "[" << SystemStats::getJUCEVersion() << "]   "
                    << String (pos.bpm, 2) << " bpm, "
                    << pos.timeSigNumerator << '/' << pos.timeSigDenominator
                    << "  -  " << timeToTimecodeString (pos.timeInSeconds)
                    << "  -  " << quarterNotePositionToBarsBeatsString (pos.ppqPosition,
                                                                        pos.timeSigNumerator,
                                                                        pos.timeSigDenominator);

        if (pos.isRecording)
            displayText << "  (recording)";
        else if (pos.isPlaying)
            displayText << "  (playing)";

        timecodeDisplayLabel.setText (displayText.toString(), dontSendNotification);
    }
}
