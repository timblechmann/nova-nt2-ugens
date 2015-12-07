HyperbolSaturation : PureUGen {
    *ar { |sig level=1|
        ^this.multiNew( 'audio', sig, level )
    }
    *kr { |sig level=1|
        ^this.multiNew( 'control', sig, level )
    }

    checkInputs { ^this.checkNInputs(1) }
}

ParabolSaturation : HyperbolSaturation {
}

PowSaturation : HyperbolSaturation {
}

TanhSaturation : HyperbolSaturation {
}

FastTanhSaturation : HyperbolSaturation {
}

FasterTanhSaturation : HyperbolSaturation {
}

NovaLinearWaveFolder : HyperbolSaturation {
}

NovaSinWaveFolder : HyperbolSaturation {
}

NovaSinWaveFolderRaw : HyperbolSaturation {
}

NovaLinearWaveWrapper : HyperbolSaturation {
}
