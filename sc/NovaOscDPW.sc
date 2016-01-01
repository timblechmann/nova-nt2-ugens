NovaSawDPW : PureUGen {
    *ar { arg freq = 440, initialPhase = 0;
        ^this.multiNew( 'audio', freq, initialPhase )
    }

    *kr { arg freq = 1, initialPhase = 0;
        ^this.multiNew( 'control', freq, initialPhase )
    }
}

NovaTriDPW : NovaSawDPW {}


NovaPulseDPW : PureUGen {
    *ar { arg freq = 440, width = 0.5, initialPhase = 0;
        ^this.multiNew( 'audio', freq, width, initialPhase )
    }

    *kr { arg freq = 1, width = 0.5, initialPhase = 0;
        ^this.multiNew( 'control', freq, width, initialPhase )
    }
}
