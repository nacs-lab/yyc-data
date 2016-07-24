# Problem and possible solutions

There are two different ways to do Raman sideband cooling AFAIK.

* CW, i.e. apply optical pumping and Raman at the same time

    The advantage is that the Raman Rabi frequency doesn't have to be
    precisely calibrated (or have long coherence time at all on the sideband).

* Pulse, i.e. drive a ~pi pulse on the sideband and do optical pumping after that

    The advantage is that the shape of the Raman pulse can be precisely
    controlled to suppress off resonance coupling on the carrier.

For sodium, the main difficulty we're having for Raman cooling is that
we have not been able to get a very clean Raman spectrum with
non-copropagating Raman beams. Whether or not there's axial component,
there always seems to be multiple closely spaced resonances near the carrier
and the resonances continuous without a clear gap into the frequency we expect
the radial sideband. This makes it hard to measure the temperature or driving
a transition on the carrier.

Because of this, for CW cooling, we don't have a good temperature measure
that is only sensitive to the axis we are cooling, combining with the strong
coupling to the cycling transition and therefore the missing of a dark state,
it is hard to measure the progress of the cooling and optimize the performance.
For pulse cooling, we can't really drive a simple pulse on the sideband that
has a high enough transfer ratio.

To help with these issues, we can make small modifications to the two approaches.
The detection is actually a common problem for both approach,
the method we can use to improve it is by doing multiple cycles of cooling and
hope to see some accumulated effect.

* For CW cooling, we can measure the progress and the fidelity of the transfer
  by discretize the process. By removing the F2 from the OP, we can avoid the
  heating due to coupling to the cycling transition and also measure the
  progress using the population accumulating in the (2,-1) state.
  This process can be repeated on the 2,-1 state by applying Raman transition
  on the (2,-1) state instead of (2,-2).

* For pulse cooling, we can improve the transfer ratio by applying a non-square
  pulse as well as sweeping the frequency of the Raman laser (LZ sweep).

# Current progress

The attempt to do pulse cooling/transfer with a LZ sweep does give better
transfer ratio although the spectrum is still hard (harder) to interpreted.
It is likely that we can improve the signal here,
now that we can afford to take more data.

The RP only OP + Raman CW method have shown good population transfer with
co-propagating Raman beams.
With non-co-propagating Raman, there seems to be some leakage out of the
(2,-1) state (most likely into other F2 states). It's unclear what the leaking
mechanism is since it doesn't happen (or much more slowly) with co-propagating
beams. Possible causes and things to try/study includes,

* Off resonance scattering from the Raman beams

    This can be better optimized/controlled with the wavemeter now.
    The leakage seems to have some detuning dependency.

* AC trap

    Might cause unexpected excitation. OP might not work with DC trap.
    The OP intensity can probably be increase to provide enough scattering
    rate in DC trap but the polarization can be messy and the power broadening
    might cause unwanted coupling to the F2 state. The F1 and F2 are merely
    176 linewidth separated and the off resonance coupling could be an issue
    with a high saturation parameter. We need a really low scattering rate
    (no more than 100kHz comparable to the Rabi rate) so this might not matter
    too much.

* DC OP

    The component when the trap is on might cause coupling with wrong
    polarization. This effect should be small with sodium and can be
    improved by increasing the B field.

* Coupling to wrong Zeeman level

    This is the main disadvantage I can think of comparing to a conventional
    CW Raman cooling. Due to the population accumulation in 2,-1, the atom
    can see of resonance Raman that is tuned for 2,-2. The 2,-2 cooling is also
    heating for 2,-1.

    The width of the OP broadened Raman resonance doesn't seem to be wide enough
    for this to be a big issue. This can also be helped by increasing B field
    and therefore the Zeeman splitting. For the first cycle, reversing the B
    field (i.e. pump to 2,2) can also help since the cooling for 2,2 is also
    cooling for 2,1.
