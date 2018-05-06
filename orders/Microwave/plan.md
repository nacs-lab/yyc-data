# Requirements

The transition we'd like to drive is from Cs (F=3, mF=3) to Cs (F=4, mF=4).
Since the binding energy of the molecular state that asymptote to (F=4, mF=4) is
298MHz, the frequency we need should be about 8.895GHz (plus Zeeman shift of up to 20MHz).
We also need to drive the atomic transition
(at least for calibration and potentially for actual state transfer in the experiment)
so we need to hit 9.193GHz (plus Zeeman shift of up to 20MHz) as well.

It's unclear how much power we actually need but a few Watts should certainly give us
an atomic resonance and there's a high chance that we can see a molecular resonance too
if we wait long enough. This also seems to be the power above
which the price and availability of amplifiers changes dramatically.

# Signal generation.

The theoretical limit what the DDS can output is 1.75GHz (half the clock) although
realistically anything higher than around 1.5GHz will be hard to obtain without sidebands.
We certainly need some frequency mixing/multiplying and here are some options.

1. Direct multiplying with a PLL

    There's chip that does exactly this from Analog Devices, e.g.
    http://www.analog.com/en/design-center/evaluation-hardware-and-software/evaluation-boards-kits/EVAL-ADF5356.html
    It's very similar to the board we use for locking the Cs laser and it can generate
    an output up to 13.6GHz at 2dBm.
    (Some version of the eval board we are using also have an VCO option but
    according to http://www.analog.com/media/en/technical-documentation/user-guides/UG-383.pdf
    the output range is 11.5GHz to 12.5GHz).

    This is the solution that requires the fewest part and is also quite cheap ($300).
    We also don't need to worry about filtering out sideband ourselves as much.

2. Mixing with a fixed source

    There are a few options for the fixed source. The worth noting ones are,

    1. Valon 5009

        This can go to 6GHz (~$250 per channel). It needs to be mixed with DDS and doubled
        to get 9GHz (doubling and then mix with DDS works too).

    2. ZSN-7800A+

        From minicircuits (https://www.minicircuits.com/pdfs/ZSN-7800A+.pdf, $500)
        This can go up to 7800MHz so we only need one mixing.

3. Direct multiplying with doublers/multipliers

    This is similar to multiplying with PLL but can allow (very non-linear) amplitude modulation
    at the same time. It removes the need for a RF switch for using a PLL + VCO which
    is surprisingly hard to find.

    Possible plan:

    DDS: 1.1 - 1.165 GHz (9dBm)
    Filter: ZX75BP-1100+ (8.5dBm) $59.95
    Double: ZX90-2-19+ (-2.5dBm) $35.95
    Filter: VBF-2275+ (-4.5dBm) $34.95
    Amp: ZX60-P103LN+ (9dBm) $69.95
    Double: ZX90-2-36-S+ (-1dBm) $36.95
    Filter: VBF-4440+ (-2dBm) $34.95
    Amp: ZX60-V63+ (18dBm) (May need attenuator or more filters) $49.95
    Double: ZX90-2-50+ (-2dBm) $41.95
    Filter: VHF-7150+ (-3dBm) $24.95
    Amp: ZX60-183A+ (25dBm) $169.95

# Pre-amp

The necessity of this depends on the source and amplifier.
There are a lot of options from minicircuits and we have a few for Cs lock so
nothing we need to buy here.

# Power Amplifier

The best one I've found is (https://www.pasternack.com/9.5-ghz-high-power-amplifier-23.5-db-gain-ip3-8.5-db-sma-pe15a5009-p.aspx) which outputs 6W. Minicircuits also has one (https://www.minicircuits.com/pdfs/ZVE-3W-183+.pdf) that's 1.4k at 3W.

Most of the amplifiers in the 3-10W range are around 2k and these two are the cheapest
among the ones of similar output power in the frequency range.

# Impedance matching with the antenna

From KRb's experience, I strongly prefer to not do this on the antenna.
Especially for near field use, it seems that the environment will strongly affect
the reflection so off-line tuning may not work at all.
It should be possible however to do this completely external to the antenna
and that, I believe, is what impedance tuners do.

The document about this in the microwave range is really scarce AFAICT.
If we go with a cut waveguide as the antenna (which is what I plan to do to start at least)
this may not be strictly needed.
I did buy two microwave impedance tuners from wester test systems
that will arrive on monday for ~$300 each although I don't really know how to understand
their performance figures. We can at least test them on the side to see how much they
can reduce the reflection and what reflection they can deal with.

# Antenna

Assuming the impedance matching issue can be done externally,
the only thing we need to consider for the antenna is the coupling to the atom.
This include the direction, polarization and distance.

The two options I come up with are cut waveguide from a 45deg angle down and
a single loop next to the chamber opposite the objective.
The waveguide wins in direction (it doesn't radiate back) tie in polarization
(both of them can realize linear polarization perpendicular to the B field) and loss
in distance.

The direction issue of the single loop antenna can potentially be fixed by placing an
reflector behind the loop (or maybe just a ground plane without a hole?) at the cost
of optical access and the distance issue of cut waveguide could be harder to fix.
It's almost certainly better to start with the waveguide and if it works we could explore
the single loop in the future (the gain shouldn't be too big though so maybe not...)
