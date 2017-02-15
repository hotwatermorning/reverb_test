//
//  Reverb.cpp
//  JuceDemoPlugin
//
//  Created by yuasa on 2017/01/23.
//
//

#include "../JuceLibraryCode/AppConfig.h"
#include "Reverb.hpp"
#include <vector>
#include <cassert>

#pragma once

#include <cmath>

//! Panpot値によって変化する音量レベルを取得する。
/*!
 "Downloadable Sounds Specification Level 2"のP.22に定義されている
 sin/cosカーブによる等パワー曲線で計算される音量レベルを返す。
 @param position Panpotの振り具合を指定する。-1.0がLeftに振り切った状態。1.0がRightに振り切った状態。0でCenterを表す。
 @return positionに指定された振り具合によって計算されるLeft側の音量レベルを返す。
 positionが-1.0(Left)の時は+3dBにあたる1.4142...が返る。
 positionが0(Center)の時は、0dBにあたる1.0が返る。
 positionが1.0(Right)の時は、-InfdBにあたる0.0が返る。
 Right側の音量レベルは、Centerを中心にしてLeftと線対称なので、positionの値に-1を掛けて振り具合を逆に指定することで取得できる。
 @note そのままの定義ではLからRまでの音量レベルを揃えるためにCenter位置でLRそれぞれの音量レベルが-3dBされてしまうので、
 それを補正するために、全体の音量を+3dBしている。
 */
inline
double GetPanVolumeBySin(double position)
{
    if(position == 1.0) {
        return 0;
    }
    return cos(M_PI * (position + 1) / 4) / cos(M_PI / 4.0);
}

//! Panpot値によって変化する音量レベルを取得する。
/*!
 "Downloadable Sounds Specification Level 1"のP.15に定義されている
 sqrtカーブによる等パワー曲線で計算される音量レベルを返す。
 Centerでの音量レベルを0dBにするため、全体の音量を+3dBしている。
 @param position Panpotの振り具合を指定する。-1.0がLeftに振り切った状態。1.0がRightに振り切った状態。0でCenterを表す。
 @return positionに指定された振り具合によって計算されるLeft側の音量レベルを返す。
 positionが-1.0(Left)の時は+3dBにあたる1.4142...が返る。
 positionが0(Center)の時は、0dBにあたる1.0が返る。
 positionが1.0(Right)の時は、-InfdBにあたる0.0が返る。
 Right側の音量レベルは、Centerを中心にしてLeftと線対称なので、positionの値に-1を掛けて振り具合を逆に指定することで取得できる。
 @note そのままの定義ではLからRまでの音量レベルを揃えるためにCenter位置でLRそれぞれの音量レベルが-3dBされてしまうので、
 それを補正するために、全体の音量を+3dBしている。
 */
inline
double GetPanVolumeBySqrt(double position)
{
    if(position == 1.0) {
        return 0;
    }
    
    return sqrt((-position + 1.0) / 2.0) / sqrt(0.5);
}

inline
double GetPanpotPositionBySin(double volume)
{
    return (4 * acos(volume * cos(M_PI / 4.0)) / M_PI) - 1;
}

inline
double linear_to_dB_abs(double linear)
{
    static double const dB_640 = 0.00000000000000000000000000000001;
    if(linear < dB_640) {
        return -640;
    } else {
        return 20.0 * log10(linear);
    }
}

//! 線形な音量値からdBへの変換
//! linearがマイナスの場合はlinearの絶対値が使用される。
//! 4.0 -> +   12dB
//! 2.0 -> +   6dB
//! 1.0 -> +/- 0dB
//! 0.5 -> -   -6dB
//! 0.25-> -   -12dB
//! 0.0 -> -   inf.dB
inline
double linear_to_dB(double linear)
{
    return linear_to_dB_abs(fabs(linear));
}

//! dBから線形な音量値への変換
//! +  12dB -> 4.0
//! +   6dB -> 2.0
//! +/- 0dB -> 1.0
//! -   6dB -> 0.5
//! -  12dB -> 0.25
//! -inf.dB-> 0.0
inline
double dB_to_linear(double dB)
{
    return pow(10.0, dB/20.0);
}

typedef float sample_t;

struct OnePole
{
    //! 正規化周波数でカットオフ周波数を設定
    /*! @pre 0 < f < 0.5
     */
    OnePole(double f)
    {
        double const k = cos(2 * M_PI * f);
        a_ = 2.0 - k - sqrt(pow(k - 2, 2) - 1);
        b_ = -a_;
    }
    
    OnePole(double a, double b)
    {
        a_ = a;
        b_ = b;
    }
    
    sample_t process(sample_t x)
    {
        x *= b_;
        x += z1_ * (-a_);
        z1_ = x;
        return x;
    }
    
private:
    sample_t a_;
    sample_t b_;
    sample_t z1_;
};

struct DelayLine
{
    DelayLine(int length)
    :   line_(length)
    ,   index_(0)
    {
    }
    
    void tick(sample_t x)
    {
        line_[index_] = x;
        index_ += 1;
        if((index_ + 1) >= line_.size()) {
            index_ = 0;
        }
    }
    
    //! get delayed sample
    sample_t get() const
    {
        return line_[index_];
    }
    
    size_t size() const
    {
        return line_.size();
    }
    
    std::vector<sample_t> line_;
    int index_;
    sample_t last_out_;
};

struct FFCF
{
    FFCF(sample_t g, int tap)
    :   g_(g)
    ,   delay_line_(tap)
    {}
    
    sample_t process(sample_t x)
    {
        sample_t const xm = delay_line_.get();
        delay_line_.tick(x);
        return x + xm * g_;
    }
    
    sample_t g_;
    DelayLine delay_line_;
};

struct FBCF
{
    FBCF(sample_t g, int tap)
    :   g_(g)
    ,   delay_line_(tap)
    {}
    
    sample_t process(sample_t x)
    {
        auto xm = delay_line_.get();
        x += xm * g_;
        delay_line_.tick(x);
        return xm;
        // ここで x + xmとすると、ccrmaの定義になる
    }
    
    sample_t g_;
    DelayLine delay_line_;
};

//! Filtered Back Feedback Comb Filter
struct FFBCF
{
    FFBCF(double cutoff, sample_t g, int tap)
    :   g_(g)
    ,   delay_line_(tap)
    ,   cutoff_(0)
    {
        SetCutoff(cutoff);
    }
    
    FFBCF(FFBCF const &rhs)
    :   delay_line_(rhs.delay_line_.size())
    {
        cutoff_ = -1;
        g_ = rhs.g_;
        delay_line_ = rhs.delay_line_;
        SetCutoff(rhs.cutoff_);
    }
    
    FFBCF & operator=(FFBCF const &rhs)
    {
        cutoff_ = -1;
        g_ = rhs.g_;
        delay_line_ = rhs.delay_line_;
        SetCutoff(rhs.cutoff_);
        
        return *this;
    }
    
    void SetCutoff(double cutoff)
    {
        if(cutoff_ != cutoff) {
            auto coeff = IIRCoefficients::makeLowPass(44100, cutoff);
            iir_.setCoefficients(coeff);
            cutoff_ = cutoff;
        }
    }
    
    sample_t process(sample_t x)
    {
        auto xm = delay_line_.get();
        x += iir_.processSingleSampleRaw(xm * g_);
        delay_line_.tick(x);
        return xm;
        // ここで x + xmとすると、ccrmaの定義になる
    }
    
    double cutoff_;
    sample_t g_;
    DelayLine delay_line_;
    IIRFilter iir_;
};

struct APF
{
    APF(sample_t g, int tap)
    :   g_(g)
    ,   delay_line_(tap)
    {}
    
    sample_t process(sample_t x)
    {
        auto xm = delay_line_.get();
        x += xm * g_;
        delay_line_.tick(x);
        return xm + x * (-g_);
    }
    
    sample_t g_;
    DelayLine delay_line_;
};

std::vector<APF> apfs = {
    { 0.7, 347 },
    { 0.7, 113 },
    { 0.7, 37 },
};

std::vector<FFBCF> ffbcfs = {
    { 6000, 0.773, 1687 },
    { 6000, 0.802, 1601 },
    { 6000, 0.753, 2053 },
    { 6000, 0.733, 2251 },
    { 6000, 0.773, 1597 },
    { 6000, 0.802, 1781 },
    { 6000, 0.753, 1771 },
    { 6000, 0.733, 2069 },
};

struct MultiTapDelay
{
    MultiTapDelay()
    {}
    
    void addDelay(sample_t g, int tap)
    {
        delays_.push_back(FFCF(g, tap));
    }
    
    sample_t process(sample_t x)
    {
        sample_t sum = 0;
        for(auto &delay: delays_) {
            sum += delay.process(x);
        }
        
        return sum / delays_.size();
    }
    
private:
    std::vector<FFCF> delays_;
    
    void CopyDelay(FFCF const &src, FFCF &dest)
    {
        dest.g_ = src.g_;
        dest.delay_line_.last_out_ = src.delay_line_.last_out_;
        if(src.delay_line_.size() > dest.delay_line_.size()) {
            auto diff = src.delay_line_.size() - dest.delay_line_.size();
            std::copy(src.delay_line_.line_.begin() + diff,
                      src.delay_line_.line_.end(),
                      dest.delay_line_.line_.begin());
            dest.delay_line_.index_ = src.delay_line_.index_ - diff;
        } else {
            std::fill(dest.delay_line_.line_.begin(),
                      dest.delay_line_.line_.end(),
                      0.0);
            
            auto diff = dest.delay_line_.size() - src.delay_line_.size();
            std::copy(src.delay_line_.line_.begin(),
                      src.delay_line_.line_.end(),
                      dest.delay_line_.line_.begin() + diff);
            dest.delay_line_.index_ = src.delay_line_.index_ - diff;
        }
    }
};

MultiTapDelay CreateEarlyReflection(float milliseconds)
{
    MultiTapDelay md;
    md.addDelay(0.7, 44100 * milliseconds / 1000.0);
    md.addDelay(0.3, 44100 * milliseconds * 1.13 / 1000.0);
    md.addDelay(0.1, 44100 * milliseconds * 1.37 / 1000.0);
//    md.addDelay(0.4, 44100 * milliseconds * 1.89 / 1000.0);
//    md.addDelay(0.32, 44100 * milliseconds * 2.39 / 1000.0);
//    md.addDelay(0.21, 44100 * milliseconds * 2.79 / 1000.0);
//    md.addDelay(0.10, 44100 * milliseconds * 3.19 / 1000.0);
    
    return md;
}

MultiTapDelay CreateEarlyReflection(float milliseconds);

struct TestReverb
{
    TestReverb(bool left)
    {
        md = CreateEarlyReflection(10);
        
        if(left) {
            apfs = {
                { 0.7, 1051 },
                { 0.7, 347 },
                { 0.7, 113 }
            };
            
            ffbcfs = {
                { 6000, 0.773, 1687 },
                { 6000, 0.802, 1601 },
                { 6000, 0.753, 2053 },
                { 6000, 0.733, 2251 },
            };
        } else {
            apfs = {
                { 0.7, 1027 },
                { 0.7, 361 },
                { 0.7, 157 }
            };
            
            ffbcfs = {
                { 6000, 0.773, 1691 },
                { 6000, 0.802, 1597 },
                { 6000, 0.753, 2009 },
                { 6000, 0.733, 2387 },
            };
        }
        
        SetErDamping(5000);
        SetPreDelayTime(20);
    }
    
    //! 周波数ごとの利得を求めて、最も周波数が大きいものに合わせてgain調整を行う
    //! z変換して周波数特性を求める
    //! 組み合わせたときの周波数特性の計算方法を調べる
    
    sample_t process(sample_t new_x, double er_gain)
    {
        // early reflection
        sample_t er = md.process(new_x);
        er = er_damping_.processSingleSampleRaw(er);

        sample_t x = 0;
        
        if(pre_delay_.size() != 0) {
            x = pre_delay_.get();
            pre_delay_.tick(new_x);
        } else {
            x = new_x;
        }
        
        for(auto &apf: apfs) {
            x = apf.process(x);
        }
        
        sample_t tmp = 0;
        
        for(int n = 0; n < ffbcfs.size(); ++n) {
            tmp += ffbcfs[n].process(x);
        }
        
        tmp /= (ffbcfs.size() + er_gain + gain_);
        
        last_sample_ = tmp + (er * er_gain);
        return last_sample_;
    }
    
    void SetTimeParam(float seconds)
    {
        auto calc_g = [seconds](auto const &filter) {
            return pow(10.0, -3 * (filter.delay_line_.size() / 44100.0) / seconds);
        };
        
        for(auto &f: ffbcfs) {
            f.g_ = calc_g(f);
        }
        
        gain_ = seconds;
    }
    
    void SetErTime(float milliseconds)
    {
        md = CreateEarlyReflection(milliseconds);
    }
    
    void SetAirAbsorption(float normalized_freq)
    {
        for(auto &f: ffbcfs) {
            f.SetCutoff(normalized_freq);
        }
    }
    
    void SetErDamping(float frequency)
    {
        er_damping_.setCoefficients(IIRCoefficients::makeLowPass(44100, frequency));
    }
    
    void SetPreDelayTime(float millisec)
    {
        pre_delay_ = DelayLine((int)std::round(millisec * 44100 / 1000.0));
    }
    
    MultiTapDelay md;
    std::vector<APF> apfs;
    std::vector<FFBCF> ffbcfs;
    IIRFilter er_damping_;
    double last_sample_ = 0;
    double gain_ = 1.0;
    
    DelayLine pre_delay_ = {441};
};

TestReverb rev_left_(true);
TestReverb rev_right_(false);

void SetTimeParam(float seconds)
{
    rev_left_.SetTimeParam(seconds);
    rev_right_.SetTimeParam(seconds);
}

void SetErTime(float milliseconds)
{
    rev_left_.SetErTime(milliseconds);
    rev_right_.SetErTime(milliseconds);
}

void SetAirAbsorption(float freq)
{
    rev_left_.SetAirAbsorption(freq);
    rev_right_.SetAirAbsorption(freq);
}

void SetErDamping(float frequency)
{
    rev_left_.SetErDamping(frequency);
    rev_right_.SetErDamping(frequency);
}

void SetPreDelayTime(float milliseconds)
{
    rev_left_.SetPreDelayTime(milliseconds);
    rev_right_.SetPreDelayTime(milliseconds);
}

double wet1 = 0.7;
double wet2 = 1.0 - wet1;

void SetStereoSpread(float amount)
{
    wet1 = (amount + 1.0) / 2.0;
    wet2 = 1.0 - wet1;
}

DelayLine fb_delay_from_left_(3720);
DelayLine fb_delay_from_right_(3163);
float fb_decay_ = 0.3;

struct DitorroReverb
{
    DelayLine pre_delay_;
    IIRFilter band_width_;
    std::vector<APF> initial_apfs_;
    struct ChReverb
    {
        APF apf1_;
        double decay_diffusion1_;
        DelayLine delay_line_;
        IIRFilter lpf_;
        double decay_;
        APF apf2_;
        double decay_diffusion2_;
        DelayLine feedback_delay_;
        double feedback_decay_;
    };
};

void ApplyReverb(AudioBuffer<float> const & input,
                 AudioBuffer<float> & output,
                 int length,
                 float dry_wet,
                 float er_gain)
{
    for(int i = 0; i < length; ++i) {
//        sample_t left = input.getReadPointer(0)[i] * (0.9 - fb_decay_) + fb_delay_from_right_.get() * fb_decay_;
//        sample_t right = input.getReadPointer(1)[i] * (0.9 - fb_decay_) + fb_delay_from_left_.get() * fb_decay_;
//        
//        left = rev_left_.process(left, dry_wet, er_gain);
//        right = rev_right_.process(right, dry_wet, er_gain);
//        
//        output.getWritePointer(0)[i] = left;
//        output.getWritePointer(1)[i] = right;
//        
//        fb_delay_from_left_.tick(left / );
//        fb_delay_from_right_.tick(right);
        
        double const dry_gain = GetPanVolumeBySin(dry_wet * 2 - 1) * cos(M_PI / 4.0);
        double const wet_gain = GetPanVolumeBySin(1 - dry_wet * 2) * cos(M_PI / 4.0);
        
        sample_t dry_left = input.getReadPointer(0)[i];
        sample_t dry_right = input.getReadPointer(1)[i];

        sample_t wet_left = rev_left_.process(dry_left, er_gain);
        sample_t wet_right = rev_right_.process(dry_right, er_gain);
        
        sample_t out_left = dry_left * dry_gain + wet_gain * (wet1 * wet_left + wet2 * wet_right);
        sample_t out_right = dry_right * dry_gain + wet_gain * (wet2 * wet_left + wet1 * wet_right);
        
        output.getWritePointer(0)[i] = out_left;
        output.getWritePointer(1)[i] = out_right;
    }
}