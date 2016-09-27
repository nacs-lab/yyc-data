# Plan to get single Cs and Na atoms simultaneously

## Current problems

1. With Na dispenser on, Cs background is too high to see live signal

    Only a problem for live signal

2. Can't distinguish two atoms since they are on the same pixel

3. Tweezer for one species can kill the loading for the other species

## Things to try

The Cs background pressure is really hard to solve since we can't control the
amount of Cs coming off the Na dispenser. It is possible to see a Cs live
signal by making the Cs MOT worse. However, doing this is unlikely going to
help anything else. If trying to get dual atom image after dropping the MOT
does not work out, we could give this a try.

All the other issues we've seen so far are due to the two tweezers overlaps too
well and the solution is to simply shift them relative to each other.
The main difficulty for solving this problem is to avoid the aberration when
we move the beam to hit different parts of the optics.
The ultimate solution to this problem is to use the AOM that we want to have
in the final setup to move the traps and also use the relay setup that is
designed to solve the aberration problem when moving the trap in the whole view.
However, that requires a full beam path reconstruction that will probably take
a few weeks and we'd like to construct and test it offline first to minimize the
down time of the experiment.

Therefore the short term solution that we'll try is to simply moving the beams
with the mirrors. Since we don't need to move it by very much, we might be
able to live with the small aberration. There are a few ways to move the beam,

1. Move the Na tweezer with 1 mirror

    This is the simplest thing to try. The issue is that we might get
    aberration. So we can,

2. Compensate the movement with a second mirror

    This should get us pretty close. There might not be space for a second
    mirror in the Na beam path so we can try instead

3. Shift the position of the beam back with a glass plate

    This also depends on if we have enough space to put in a glass plate
    but might be easier to implement since it doesn't need changing the
    beam path geometry. We can potentially temporarily remove some
    optics from the beam path to make space for the glass plate or we might
    even,

4. Use the PBS that's already in the beam path to shift the beam.

    We obviously can't change the angle of the PBS too much but we could try.

5. We can also shift the Cs tweezer, which has two mirrors
