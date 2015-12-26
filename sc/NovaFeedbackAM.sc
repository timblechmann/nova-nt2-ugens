NovaFeedbackAM : PureUGen {
        *ar { arg sig, fb;
                ^this.multiNew( 'audio', sig, fb )
        }

        *kr { arg sig, fb;
                ^this.multiNew( 'control', sig, fb )
        }
}

NovaFeedbackAM2 : MultiOutUGen {
        *ar { arg left, right, fb;
                ^this.multiNew( 'audio', left, right, fb )
        }

        *kr { arg left, right, fb;
                ^this.multiNew( 'control', left, right, fb )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..2).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(2) }
}

NovaFeedbackAM2_2 : MultiOutUGen {
        *ar { arg left, right, fb0, fb1;
                ^this.multiNew( 'audio', left, right, fb0, fb1 )
        }

        *kr { arg left, right, fb0, fb1;
                ^this.multiNew( 'control', left, right, fb0, fb1 )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..2).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(2) }
}

NovaFeedbackAM4 : MultiOutUGen {
        *ar { arg s0, s1, s2, s3, fb;
                ^this.multiNew( 'audio', s0, s1, s2, s3, fb )
        }

        *kr { arg s0, s1, s2, s3, fb0, fb;
                ^this.multiNew( 'control', s0, s1, s2, s3, fb)
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..3).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(4) }
}

NovaFeedbackAM4_4 : MultiOutUGen {
        *ar { arg s0, s1, s2, s3, fb0, fb1, fb2, fb3;
                ^this.multiNew( 'audio', s0, s1, s2, s3, fb0, fb1, fb2, fb3 )
        }

        *kr { arg s0, s1, s2, s3, fb0, fb1, fb2, fb3;
                ^this.multiNew( 'control', s0, s1, s2, s3, fb0, fb1, fb2, fb3 )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..3).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(4) }
}

NovaFeedbackAM8_8 : MultiOutUGen {
        *ar { arg s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7;
                ^this.multiNew( 'audio', s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7 )
        }

        *kr { arg s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7;
                ^this.multiNew( 'control', s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7 )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..7).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(8) }
}
