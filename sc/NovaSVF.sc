//////////////////////
//
// svf based filters
//

NovaLowPassSVF : Filter {
    *ar { |sig, freq = 880, res = 0|
        ^this.multiNew('audio', sig, freq, res)
    }
}

NovaHighPassSVF   : NovaLowPassSVF {}
NovaBandPassSVF   : NovaLowPassSVF {}
NovaBandRejectSVF : NovaLowPassSVF {}
NovaPeakSVF       : NovaLowPassSVF {}

NovaLowShelfSVF : Filter {
    *ar { |sig, freq = 880, amp = 1, res = 0|
        ^this.multiNew('audio', sig, freq, amp, res)
    }
}

NovaHighShelfSVF  : NovaLowShelfSVF {}
NovaEqSVF         : NovaLowShelfSVF {}



/////////////////////////
// 2-channel, scalar args

NovaLowPassSVF2 : PureMultiOutUGen {
    *ar { |sig1, sig2, freq = 880, res = 0|
        ^this.multiNew('audio', sig1, sig2, freq, res)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}

NovaHighPassSVF2   : NovaLowPassSVF2 {}
NovaBandPassSVF2   : NovaLowPassSVF2 {}
NovaBandRejectSVF2 : NovaLowPassSVF2 {}
NovaPeakSVF2       : NovaLowPassSVF2 {}

NovaLowShelfSVF2 : PureMultiOutUGen {
    *ar { |sig1, sig2, freq = 880, amp = 1, res = 0|
        ^this.multiNew('audio', sig1, sig2, freq, amp, res)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }

}

NovaHighShelfSVF2  : NovaLowShelfSVF2 {}
NovaEqSVF2         : NovaLowShelfSVF2 {}


///////////////////////////
// 2-channel, separate args

NovaLowPassSVF2_2 : PureMultiOutUGen {
    *ar { |sig1, sig2, freq1 = 880, freq2 = 880, res1 = 0, res2 = 0|
        ^this.multiNew('audio', sig1, sig2, freq1, freq1, res1, res2)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}

NovaHighPassSVF2_2   : NovaLowPassSVF2_2 {}
NovaBandPassSVF2_2   : NovaLowPassSVF2_2 {}
NovaBandRejectSVF2_2 : NovaLowPassSVF2_2 {}
NovaPeakSVF2_2       : NovaLowPassSVF2_2 {}

NovaLowShelfSVF2_2 : PureMultiOutUGen {
    *ar { |sig1, sig2, freq1 = 880, freq2 = 880, amp1 = 1, amp2 = 1, res1 = 0, res2 = 0|
        ^this.multiNew('audio', sig1, sig2, freq1, freq1, amp1, amp2, res1, res2)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }

}

NovaHighShelfSVF2_2  : NovaLowShelfSVF2_2 {}
NovaEqSVF2_2         : NovaLowShelfSVF2_2 {}


/////////////////////////
// 4-channel, scalar args

NovaLowPassSVF4 : PureMultiOutUGen {
    *ar { |sig1, sig2, sig3, sig4, freq = 880, res = 0|
        ^this.multiNew('audio', sig1, sig2, sig3, sig4, freq, res)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..3).collect({|i| OutputProxy(rate, this, i) });
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}

NovaHighPassSVF4   : NovaLowPassSVF4 {}
NovaBandPassSVF4   : NovaLowPassSVF4 {}
NovaBandRejectSVF4 : NovaLowPassSVF4 {}
NovaPeakSVF4       : NovaLowPassSVF4 {}

NovaLowShelfSVF4 : PureMultiOutUGen {
    *ar { |sig1, sig2, sig3, sig4, freq = 880, amp = 1, res = 0|
        ^this.multiNew('audio', sig1, sig2, sig3, freq, amp, res)
    }
}

NovaHighShelfSVF4  : NovaLowShelfSVF4 {}
NovaEqSVF4         : NovaLowShelfSVF4 {}
