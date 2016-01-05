NovaTiltFilter : PureUGen {
    *ar { arg sig, cutoff = 440, gain = 0;
        ^this.multiNew( 'audio', sig, cutoff, gain )
    }

    *kr { arg sig, cutoff = 10, gain = 0;
        ^this.multiNew( 'control', sig, cutoff, gain )
    }
}

NovaTiltFilter2 : PureMultiOutUGen {
    *ar { arg left, right, cutoff = 440, gain = 0;
        ^this.multiNew( 'audio', left, right, cutoff, gain )
    }

    *kr { arg left, right, cutoff = 10, gain = 0;
        ^this.multiNew( 'control', left, right, cutoff, gain )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..1).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}

NovaTiltFilter4 : PureMultiOutUGen {
    *ar { arg s0, s1, s2, s3, cutoff = 440, gain = 0;
        ^this.multiNew( 'audio', s0, s1, s2, s3, cutoff, gain )
    }

    *kr { arg s0, s1, s2, s3, cutoff = 10, gain = 0;
        ^this.multiNew( 'control', s0, s1, s2, s3, cutoff, gain )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..3).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}

NovaTiltFilter2_2 : NovaTiltFilter2 {
    *ar { arg left, right, cutoffLeft = 440, cutoffRight = 440, gainLeft = 0, gainRight = 0;
        ^this.multiNew( 'audio', left, right, cutoffLeft, cutoffRight, gainLeft, gainRight )
    }

    *kr { arg left, right, cutoffLeft = 440, cutoffRight = 440, gainLeft = 0, gainRight = 0;
        ^this.multiNew( 'control', left, right, cutoffLeft, cutoffRight, gainLeft, gainRight )
    }
}
