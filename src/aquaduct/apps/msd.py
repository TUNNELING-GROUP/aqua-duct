#coding:utf-8 -*-

# THIS IS COPY OF MDAnalysis CODE, PASTED HERE FOR TESTING PURPOSES ONLY

import numpy as np
from MDAnalysis.lib.log import ProgressMeter


class MeanSquareDisplacementAQ(object):
    r"""Mean square displacement analysis

    Function to evaluate the Mean Square Displacement (MSD_). The MSD gives the
    average distance that particles travels. The MSD is given by:

    .. math::
        \langle\Delta r(t)^2\rangle = 2nDt

    where :math:`r(t)` is the position of particle in time :math:`t`,
    :math:`\Delta r(t)` is the displacement after time lag :math:`t`,
    :math:`n` is the dimensionality, in this case :math:`n=3`,
    :math:`D` is the diffusion coefficient and :math:`t` is the time.

    .. _MSD: http://en.wikipedia.org/wiki/Mean_squared_displacement


    Parameters
    ----------
    universe : Universe
      Universe object
    selection : str
      Selection string for water [‘byres name OH2’].
    t0 : int
      frame  where analysis begins
    tf : int
      frame where analysis ends
    dtmax : int
      Maximum dt size, `dtmax` < `tf` or it will crash.


    .. versionadded:: 0.11.0
    """

    def __init__(self, universe, selection, t0, tf, dtmax, nproc=1):
        self.universe = universe
        self.selection = selection
        self.t0 = t0
        self.tf = tf
        self.dtmax = dtmax
        self.nproc = nproc
        self.timeseries = None

    def _repeatedIndex(self, selection, dt, totalFrames):
        """
        Indicate the comparation between all the t+dt.
        The results is a list of list with all the repeated index per frame
        (or time).

        - Ex: dt=1, so compare frames (1,2),(2,3),(3,4)...
        - Ex: dt=2, so compare frames (1,3),(3,5),(5,7)...
        - Ex: dt=3, so compare frames (1,4),(4,7),(7,10)...
        """
        rep = []
        for i in range(int(totalFrames-dt)):
            rep.append(self._sameMolecTandDT(
                selection, i, i + dt))
#        for i in range(int(round((totalFrames - 1) / float(dt)))):
#            if (dt * i + dt < totalFrames):
#                rep.append(self._sameMolecTandDT(
#                    selection, dt * i, (dt * i) + dt))
        return rep

    def _getOneDeltaPoint(self, universe, repInd, i, t0, dt):
        """
        Gives one point to calculate the mean and gets one point of the plot
        C_vect vs t.

        - Ex: t0=1 and dt=1 so calculate the t0-dt=1-2 interval.
        - Ex: t0=5 and dt=3 so calcultate the t0-dt=5-8 interva

        i = come from getMeanOnePoint (named j) (int)
        """
        valO = 0
        n = 0
        for j in range(len(repInd[i]) // 3):
            begj = 3 * j
            universe.trajectory[t0]
            # Plus zero is to avoid 0to be equal to 0tp
            Ot0 = repInd[i][begj].position + 0

            universe.trajectory[t0 + dt]
            # Plus zero is to avoid 0to be equal to 0tp
            Otp = repInd[i][begj].position + 0

            # position oxygen
            OVector = Ot0 - Otp
            # here it is the difference with
            # waterdynamics.WaterOrientationalRelaxation
            valO += np.dot(OVector, OVector)
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return valO/n if n > 0 else 0

    def _getMeanOnePoint(self, universe, selection1, selection_str, dt,
                         totalFrames):
        """
        This function gets one point of the plot C_vec vs t. It's uses the
        _getOneDeltaPoint() function to calculate the average.

        """
        repInd = self._repeatedIndex(selection1, dt, totalFrames)
        sumsdt = 0
        n = 0.0
        sumDeltaO = 0.0
        valOList = []

        for j in range(totalFrames // dt - 1):
            a = self._getOneDeltaPoint(universe, repInd, j, sumsdt, dt)
            sumDeltaO += a
            valOList.append(a)
            sumsdt += dt
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return sumDeltaO/n if n > 0 else 0

    def _sameMolecTandDT(self, selection, t0d, tf):
        """
        Compare the molecules in the t0d selection and the t0d+dt selection and
        select only the particles that are repeated in both frame. This is to
        consider only the molecules that remains in the selection after the dt
        time has elapsed. The result is a list with the indexs of the atoms.
        """
        a = set(selection[t0d])
        b = set(selection[tf])
        sort = sorted(list(a.intersection(b)))
        return sort

    def _selection_serial(self, universe, selection_str):
        selection = []
        pm = ProgressMeter(universe.trajectory.n_frames,
                           interval=10, verbose=True)
        l = 0
        for ts in universe.trajectory:
            selection.append(universe.select_atoms(selection_str[l]))
            pm.echo(ts.frame)
            l += 1
        return selection

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries"""

        # All the selection to an array, this way is faster than selecting
        # later.
        if self.nproc == 1:
            selection_out = self._selection_serial(
                self.universe, self.selection)
        else:
            # parallel not yet implemented
            # selection = selection_parallel(universe, selection_str, nproc)
            selection_out = self._selection_serial(
                self.universe, self.selection)
        self.timeseries = []

        pm = ProgressMeter(self.dtmax + 1,
                           interval=1, verbose=True)
        for dt in list(range(1, self.dtmax + 1)):
            output = self._getMeanOnePoint(
                self.universe, selection_out, self.selection, dt, self.tf)
            self.timeseries.append(output)
            pm.echo(dt)


if __name__ == "__main__":

    import MDAnalysis

    universe = MDAnalysis.Universe('topology','trajectory')

    with open('selection_scope.txt') as f:
        selection_AQ = f.readlines()

    MSD_analysis = MeanSquareDisplacementAQ(universe, selection_AQ, 0, 10000, 100)
    MSD_analysis.run()
