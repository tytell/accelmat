import numpy as np
from scipy import interpolate, signal
import quaternion
import h5py
from copy import copy
import re

import matplotlib.pyplot as plt

class IMU(object):
    def __init__(self):
        self.t = []
        self.gyro = []
        self.sampfreq = None

        self.Qbias = None
        self.Qgyro = None
        self.Qacc = None
        self.Qdyn = None

        self.gN = np.array([0, 0, 1])

    def set_gyro(self, t, gyro, resamplefreq=None):
        if resamplefreq is not None:
            t, gyro = self._resample_gyro(t, gyro, resamplefreq)
            sampfreq = resamplefreq
        else:
            t, gyro = t, gyro
            sampfreq = np.mean(np.diff(t))

        self.sampfreq = sampfreq
        self.t = t
        self.gyro = self._filter_gyro(t, gyro)

    def load(self, filename, t_dataset='/data/t', gyro_dataset='/data/Gyro', accel_dataset='/data/Accel',
             resamplefreq=None, t_units='ms'):
        with h5py.File(filename, 'r') as h5file:
            acc = np.array(h5file[accel_dataset])
            gyro = np.array(h5file[gyro_dataset])
            t0 = np.array(h5file[t_dataset])

            if t_units == 'ms' or t_units == 'msec':
                t0 /= 1000.0
            elif t_units == 's' or t_units == 'sec':
                pass
            else:
                raise ValueError('Unrecognized unit for time: {}'.format(t_units))

            if resamplefreq is not None:
                if isinstance(resamplefreq, basestring):
                    if resamplefreq == 'mean':
                        resamplefreq = np.around(1.0 / np.mean(np.diff(t0)), decimals=-1)
                    else:
                        raise ValueError('Unrecognized resampling option: {}'.format(resamplefreq))

                self.t_irr = t0
                self.gyro_irr = gyro
                self.acc_irr = acc

                t = np.arange(t0[0], t0[-1], 1.0/resamplefreq)

                gyro = interpolate.interp1d(t0, gyro, axis=0)(t)
                acc = interpolate.interp1d(t0, acc, axis=0)(t)
                sampfreq = resamplefreq
            else:
                t = t0
                sampfreq = 1.0/np.mean(np.diff(t0))

            self.t0 = t
            self.t = t
            self.acc0 = acc
            self.acc = self.acc0
            self.gyro0 = gyro
            self.gyro = self.gyro0
            self.sampfreq = sampfreq

    def filter(self, order=None, gyro_cutoff=None, acc_cutoff=None, nsamp=None, method='butter'):
        if method == 'butter':
            if gyro_cutoff is not None:
                if len(gyro_cutoff) == 1:
                    gyro_hi = gyro_cutoff[0]
                    gyro = self.gyro0
                elif len(gyro_cutoff) == 2:
                    gyro_hi = gyro_cutoff[1]

                    gyrolo = self.get_low_baseline(self.t0, self.gyro0, gyro_cutoff[0])
                    gyro = self.gyro0 - gyrolo
                else:
                    raise ValueError('Unrecognized frequency range {}'.format(acc_cutoff))

                gyro_sos = signal.butter(order, gyro_hi/(self.sampfreq/2.0), "lowpass", output='sos')

                self.gyro = signal.sosfiltfilt(gyro_sos, gyro, axis=0)
            else:
                self.gyro = self.gyro0

            if acc_cutoff is not None:
                if len(acc_cutoff) == 1:
                    acc_hi = acc_cutoff[0]
                    acc = self.acc0
                elif len(acc_cutoff) == 2:
                    acc_hi = acc_cutoff[1]

                    acclo = self.get_low_baseline(self.t0, self.acc0, acc_cutoff[0])
                    acc = self.acc0 - acclo
                else:
                    raise ValueError('Unrecognized frequency range {}'.format(acc_cutoff))

                acc_sos = signal.butter(order, acc_hi / (self.sampfreq / 2.0), 'lowpass', output='sos')

                self.acc = signal.sosfiltfilt(acc_sos, acc, axis=0)
            else:
                self.acc = self.acc0

        elif method == 'running':
            self.gyro = np.zeros_like(self.gyro0)
            self.accel = np.zeros_like(self.acc0)

            for i in range(3):
                self.gyro[:, i] = np.convolve(self.gyro0[:, i], np.ones((nsamp,))/nsamp, mode='same')
                self.accel[:, i] = np.convolve(self.acc0[:, i], np.ones((nsamp,))/nsamp, mode='same')

    def get_orientation(self, method='madgwick', initwindow=0.5, beta=2.86, lCa=(0.0, -0.3, -0.7), Ca=None):
        if method.lower() == 'ekf':
            dt = np.mean(np.diff(self.t))
            if Ca is None:
                Ca = np.power(10, np.array(lCa)) / dt
            else:
                Ca = np.array(Ca) / dt
            self._get_orientation_ekf(Ca=Ca)

            inertialbasis = []

            for rpy in self.orient_sensor:
                QT = self._getQT(rpy)

                inertialbasis.append(QT.dot(np.eye(3)))

            inertialbasis = np.array(inertialbasis)

            self.worldbasis = np.matmul(self.chip2world_rot, inertialbasis)

            self.accdyn_world = np.matmul(self.chip2world_rot, self.accdyn_sensor.T).T
            self.accdyn = self.accdyn_world

        elif method.lower() in ['madgwick', 'integrate_gyro']:
            if method.lower() == 'integrate_gyro':
                beta = 0.0

            self._get_orientation_madgwick(initwindow=initwindow, beta=beta)

            qorient_world = self.qchip2world.conj() * self.qorient * self.qchip2world
            self.qorient_world = qorient_world

            self.orient_world = np.array([quaternion.as_euler_angles(q1) for q1 in qorient_world])
            self.orient = self.orient_world

            # make accdyn into a quaternion with zero real part
            qacc = np.zeros((self.accdyn_sensor.shape[0], 4))
            qacc[:, 1:] = self.accdyn_sensor

            # rotate accdyn into the world coordinate system
            qaccdyn_world = self.qchip2world.conj() * quaternion.as_quat_array(qacc) * self.qchip2world
            self.accdyn_world = np.array([q.components[1:] for q in qaccdyn_world])
            self.accdyn = self.accdyn_world


        return self.orient_world

    def calibrate(self, filename):
        with h5py.File(filename, 'r') as h5calib:
            gyro = np.array(h5calib['/data/Gyro'])
            # convert file data from deg/sec to rad/sec
            gyro = np.deg2rad(gyro)

            self.bias_gyro = np.mean(gyro, axis=0)

            accel = 9.81 * np.array(h5calib['/data/Accel'])

            # get noise covariances
            self.Qgyro = np.cov(gyro, rowvar=False)
            self.Qacc = np.cov(accel, rowvar=False)

            # bias noise covariance (assuming low drift)
            self.Qbias = 1e-10 * self.Qacc

    def get_inertial_coords(self, filename, method='mean accel', g=None):
        with h5py.File(filename, 'r') as h5inertial:
            accel = 9.81 * np.array(h5inertial['/data/Accel'])

            if g is not None:
                self.gN = g
            elif method == 'mean accel':
                self.gN = np.mean(accel, axis=0)

    def get_world_coordinates(self, filename, axes=['z'], times=None, averagedur=0.1):
        axinddict = {'x': 0, 'y': 1, 'z': 2}

        with h5py.File(filename, 'r') as h5calib:
            acc = np.array(h5calib['/data/Accel'])
            t = np.array(h5calib['/data/t'])

        gax = np.eye(3)

        if len(axes) == 1:
            times = [0]
            averagedur = 4*t[-1]

        d2 = averagedur/2

        axord = []
        for axis1, time1 in zip(axes, times):
            ist = np.logical_and(t >= time1-d2, t <= time1+d2)
            ax = np.mean(acc[ist, :], axis=0)

            m = re.match('([+-]?)([XYZxyz])', axis1)
            if m is None:
                raise ValueError('Unrecognized axis specification axis')

            if m.group(1) == '-':
                axsign = -1.0
            else:
                axsign = 1.0

            axind = axinddict[m.group(2).lower()]

            gax[:, axind] = ax * axsign
            axord.append(axind)

        axord = np.concatenate((axord, np.setdiff1d(range(3), axord)))
        axrevord = np.argsort(np.arange(3)[axord])

        basis = self._gramschmidt(gax[:, axord])
        basis = basis[:, axrevord]

        # check for right handed-ness
        # the Z axis should be equal to the cross product of the X and Y axes
        # because of small numerical issues, it's sometimes not exactly equal,
        # so we check that they're in the same direction
        assert(np.dot(np.cross(basis[:, 0], basis[:, 1]), basis[:, 2]) > 0.9)

        self.chip2world_rot = basis
        self.qchip2world = quaternion.from_rotation_matrix(basis)
        self.qworld2chip = quaternion.from_rotation_matrix(basis.T)

    def _gramschmidt(self, U):
        k = U.shape[1]
        V = copy(U)

        for i in range(k):
            V[:, i] /= np.linalg.norm(V[:, i])

            for j in range(i+1, k):
                proj = np.dot(V[:, i], V[:, j]) * V[:, i]
                V[:, j] -= proj

        return V

    def _get_orientation_ekf(self, Ca):
        """Extended Kalman filter for sensor fusion
        x is the state: [theta, bias, adyn]^T
        where theta are the Euler """
        if self.Qbias is None or self.Qgyro is None or self.Qacc is None:
            raise ValueError('Noise covariance is not yet estimated')
        if self.gN is None:
            raise ValueError('Inertial reference frame is not yet set')

        self.Qdyn = np.diag(Ca).dot(self.Qacc)

        xkm1 = np.zeros((9,))
        dt = np.diff(self.t)
        dt = np.insert(dt, 0, dt[0:0])
        dt1 = dt[0]

        # make 9 x 9 matrix
        Pkm1 = self._stack_matrices([
            [(self.Qgyro + self.Qbias)*dt1**2,  -self.Qbias*dt1,    np.zeros((3, 3))],
            [-self.Qbias*dt1,                   self.Qbias,         np.zeros((3, 3))],
            [np.zeros((3, 3)),                  np.zeros((3, 3)),   self.Qdyn]])

        N = self.gyro.shape[0]

        gyro = np.deg2rad(self.gyro) - self.bias_gyro
        acc = 9.81*self.acc

        eulerEKF = []
        aD = []
        err = []
        PkEKF = []
        xkEKF = []
        Rk = self.Qacc

        for dt1, omegak, accel in zip(dt, gyro, acc):
            Fk, xkM, Bk = self._system_dynamics(xkm1, omegak, dt1, Ca)

            Qk = self._stack_matrices([
                [Bk.dot(self.Qgyro + self.Qbias).dot(Bk.T)*dt1**2,    -Bk.dot(self.Qbias)*dt1,  np.zeros((3,3))],
                [-self.Qbias.dot(Bk.T)*dt1,                           self.Qbias,               np.zeros((3,3))],
                [np.zeros((3,3)),                                     np.zeros((3,3)),          self.Qdyn]])

            PkM = Fk.dot(Pkm1).dot(Fk.T) + Qk
            hk, Jh = self._observation_dynamics(xkM, self.gN)

            Hk = Jh

            Sk = Hk.dot(PkM).dot(Hk.T) + Rk
            Kk = PkM.dot(Hk.T).dot(np.linalg.pinv(Sk))
            xk = xkM + Kk.dot(accel - hk)
            Pk = (np.eye(9) - Kk.dot(Hk)).dot(PkM)

            QT = self._getQT(xk[:3])

            eulerEKF.append(xk[:3])
            aD.append(QT.T.dot(xk[6:]))
            err.append(accel - ((QT.dot(self.gN) + xk[6:])))
            PkEKF.append(Pk)
            xkEKF.append(xk)

            Pkm1 = Pk
            xkm1 = xk

        self.orient_sensor = np.pad(np.array(eulerEKF), ((1, 0), (0, 0)), mode='edge')
        self.accdyn_sensor = np.pad(np.array(aD), ((1, 0), (0, 0)), mode='edge')

        qorient = []
        for o1 in self.orient_sensor:
            qorient.append(quaternion.from_euler_angles(*o1))
        self.qorient = np.array(qorient)

    def _system_dynamics(self, xk, omegak, dt, Ca):
        phi, theta, psi = xk[:3]
        biask = xk[3:6]

        sPh = np.sin(phi)
        cPh = np.cos(phi)
        tTh = np.tan(theta)
        scTh = 1 / np.cos(theta)

        Bk = np.array([[1,  sPh*tTh,    cPh*tTh],
                       [0,  cPh,        -sPh],
                       [0,  sPh*scTh,   cPh*scTh]])

        # partial diffs
        Bk_phi = np.array([[0,  cPh*tTh,    -sPh*tTh],
                           [0,  -sPh,       -cPh],
                           [0,  cPh*scTh,   -sPh*scTh]])

        Bk_theta = np.array([[0,    sPh*scTh**2,    cPh*scTh**2],
                             [0,    0,              0],
                             [0,    sPh*scTh*tTh,   cPh*scTh*tTh]])

        Bk_psi = np.zeros((3, 3))

        unbiased_omegak = omegak - biask
        unbiased_omegak = unbiased_omegak[:, np.newaxis]

        Bkomega = np.hstack((np.dot(Bk_phi, unbiased_omegak),
                             np.dot(Bk_theta, unbiased_omegak),
                             np.dot(Bk_psi, unbiased_omegak)))

        Fk = self._stack_matrices([
            [np.eye(3) + Bkomega*dt, -Bk*dt, np.zeros((3,3))],
            [np.zeros((3, 3)), np.eye(3), np.zeros((3, 3))],
            [np.zeros((3, 3)), np.zeros((3, 3)), np.eye(3)]])

        xkp1 = np.hstack((xk[:3] + np.dot(Bk, unbiased_omegak).squeeze()*dt, xk[3:6], xk[6:] + Ca*dt))

        return Fk, xkp1, Bk

    def _observation_dynamics(self, xk, gN):
        """gN = gravity in inertial coordinate system (3x1)"""
        phi, theta, psi = xk[:3]

        # rotation matrices
        Rz_yaw =    np.array([[np.cos(psi),     np.sin(psi),    0],
                             [-np.sin(psi),     np.cos(psi),    0],
                             [0,                0,              1]])
        Ry_pitch = np.array([[np.cos(theta),    0,              -np.sin(theta)],
                             [0,                1,              0],
                             [np.sin(theta),    0,              np.cos(theta)]])
        Rx_roll =  np.array([[1,                0,              0],
                             [0,                np.cos(phi),    np.sin(phi)],
                             [0,                -np.sin(phi),   np.cos(phi)]])

        # rates/derivatives
        Rz_yaw_rate =   np.array([[-np.sin(psi),    np.cos(psi),    0],
                                  [-np.cos(psi),    -np.sin(psi),   0],
                                  [0,               0,              0]])
        Ry_pitch_rate = np.array([[-np.sin(theta),  0,              -np.cos(theta)],
                                  [0,               0,              0],
                                  [np.cos(theta),   0,              -np.sin(theta)]])
        Rx_roll_rate =  np.array([[0,               0,              0],
                                  [0,               -np.sin(phi),   np.cos(phi)],
                                  [0,               -np.cos(phi),   -np.sin(phi)]])

        QT       = Rx_roll.dot(Ry_pitch).dot(Rz_yaw)
        QT_roll  = Rx_roll_rate.dot(Ry_pitch).dot(Rz_yaw)
        QT_pitch = Rx_roll.dot(Ry_pitch_rate).dot(Rz_yaw)
        QT_yaw   = Rx_roll.dot(Ry_pitch).dot(Rz_yaw_rate)

        Jh = np.vstack((QT_roll.dot(gN), QT_pitch.dot(gN), QT_yaw.dot(gN), np.zeros((3, 3)), np.eye(3))).T
        hk = QT.dot(gN) + xk[6:]

        return hk, Jh

    def _getQT(self, x):
        phi, theta, psi = x

        Rz_yaw =    np.array([[np.cos(psi),     np.sin(psi),    0],
                             [-np.sin(psi),     np.cos(psi),    0],
                             [0,                0,              1]])
        Ry_pitch = np.array([[np.cos(theta),    0,              -np.sin(theta)],
                             [0,                1,              0],
                             [np.sin(theta),    0,              np.cos(theta)]])
        Rx_roll =  np.array([[1,                0,              0],
                             [0,                np.cos(phi),    np.sin(phi)],
                             [0,                -np.sin(phi),   np.cos(phi)]])
        QT = Rx_roll.dot(Ry_pitch).dot(Rz_yaw)

        return QT

    def _stack_matrices(self, M):
        m = []
        for row in M:
            m.append(np.hstack(tuple(row)))
        return np.vstack(tuple(m))

    def _get_orientation_madgwick(self, initwindow=0.5, beta=2.86):
        gyrorad = self.gyro * np.pi/180.0
        betarad = beta * np.pi/180.0

        qorient = np.zeros_like(self.t, dtype=np.quaternion)
        qgyro = np.zeros_like(self.t, dtype=np.quaternion)
        gvec = np.zeros_like(self.gyro)

        dt = 1.0 / self.sampfreq

        isfirst = self.t <= self.t[0] + initwindow
        qorient[0] = self.orientation_from_accel(np.mean(self.acc[isfirst, :], axis=0))

        for i, gyro1 in enumerate(gyrorad[1:, :], start=1):
            qprev = qorient[i-1]

            acc1 = self.acc[i, :]
            acc1 = -acc1 / np.linalg.norm(acc1)

            # quaternion angular change from the gryo
            qdotgyro = 0.5 * (qprev * np.quaternion(0, *gyro1))
            qgyro[i] = qprev + qdotgyro * dt

            if beta > 0:
                # gradient descent algorithm corrective step
                qp = qprev.components
                F = np.array([2*(qp[1]*qp[3] - qp[0]*qp[2]) - acc1[0],
                              2*(qp[0]*qp[1] + qp[2]*qp[3]) - acc1[1],
                              2*(0.5 - qp[1]**2 - qp[2]**2) - acc1[2]])
                J = np.array([[-2*qp[2], 2*qp[3], -2*qp[0], 2*qp[1]],
                               [2*qp[1], 2*qp[0], 2*qp[3], 2*qp[2]],
                               [0, -4*qp[1], -4*qp[2], 0]])

                step = np.dot(J.T, F)
                step = step / np.linalg.norm(step)
                step = np.quaternion(*step)

                qdot = qdotgyro - betarad * step
            else:
                qdot = qdotgyro

            qorient[i] = qprev + qdot * dt

        # get the gravity vector
        # gravity is -Z
        gvec = [(q.conj() * np.quaternion(0, 0, 0, -1) * q).components[1:] for q in qorient]

        self.qorient = qorient
        self.orient_sensor = quaternion.as_euler_angles(qorient)
        self.gvec = gvec
        self.accdyn_sensor = self.acc - gvec

    def integrate_gyro(self, initwindow=0.5):
        # convert gyro to rad/sec
        gyrorad = self.gyro * np.pi/180.0

        qorient = np.zeros_like(self.t, dtype=np.quaternion)
        gvec = np.zeros_like(self.gyro)

        dt = 1.0 / self.sampfreq

        isfirst = self.t <= self.t[0] + initwindow
        qorient[0] = self.orientation_from_accel(np.mean(self.acc[isfirst, :], axis=0))

        for i, gyro1 in enumerate(gyrorad[1:, :], start=1):
            qprev = qorient[i-1]

            # quaternion angular change from the gyro
            qdotgyro = 0.5 * (qprev * np.quaternion(0, *gyro1))
            qorient[i] = qprev + qdotgyro * dt

            g1 = qorient[i].conj() * np.quaternion(0, 0, 0, -1) * qorient[i]
            gvec[i, :] = g1.components[1:]

        self.qorient = qorient
        self.orient_sensor = quaternion.as_euler_angles(qorient)
        self.gvec = gvec
        self.accdyn_sensor = self.acc - gvec

    def orientation_from_accel(self, acc):
        """Get a quaternion orientation from an accelerometer reading, assuming that the accelerometer correctly
        measures the gravitational acceleration"""
        acc1 = -acc / np.linalg.norm(acc)

        ax, ay, az = acc1

        AXZ = ax * np.sqrt(1 - az)
        AXY = np.sqrt(ax**2 + ay**2)
        q0 = np.quaternion(0, AXZ / (np.sqrt(2)*AXY), ay*AXZ / (np.sqrt(2) * ax * AXY), ax*AXY / (np.sqrt(2)*AXZ))

        return q0

    def get_low_baseline(self, t, y, cutoff, padmode='mean', stat_length=None):
        dur = 1.0 / cutoff
        dt = t[1] - t[0]
        n = int(dur // dt)+1
        N = len(t)

        nblocks = int(np.ceil(float(N)/n))
        print "dt={}, dur={}, N={}, n={}, nblocks={}".format(dt, dur,N,n,nblocks)

        pad = n*nblocks - N

        pad1 = int(pad // 2)
        pad2 = int(pad - pad1)

        y = np.pad(y, ((pad1, pad2), (0, 0)), mode=padmode, stat_length=stat_length)

        yblock = np.split(y, nblocks, axis=0)

        ymn = np.mean(yblock, axis=1)
        ymn = np.pad(ymn, ((1, 1), (0, 0)), mode='edge')
        ctrind = np.arange(nblocks)*n + n/2
        ctrind = np.concatenate(([0], ctrind, [N-1]))
        ind = np.arange(N)

        print "nblocks={}, ymn.shape={}".format(nblocks, ymn.shape)
        ylo = interpolate.interp1d(ctrind, ymn, kind='cubic', axis=0)(ind)
        return ylo

    def _filter_gyro(self, t, gyro):
        if len(self.freqrange) == 1:
            btype = 'lowpass'
        elif len(self.freqrange) == 2:
            btype = 'bandpass'
        else:
            raise ValueError('Unrecognized frequency range {}'.format(self.freqrange))

        sos = signal.butter(self.filterorder, self.freqrange/(self.sampfreq/2.0), btype, output='sos')

        gyros = signal.sosfiltfilt(sos, gyro, axis=0)
        return gyros

    def _resample_gyro(self, t0, gyro0, resamplefreq):
        t = np.arange(t0[0], t0[-1], 1.0/resamplefreq)

        intp = interpolate.interp1d(t0, gyro0, axis=0)
        gyro = intp(t)

        return t, gyro


def main():
    filename = '/Users/etytel01/Documents/Acceleration/rawdata/two_imu_data/data_b_2.hdf5'
    calibfilename = '/Users/etytel01/Documents/Acceleration/rawdata/two_imu_data/data_b_1.hdf5'
    # encoderfilename = '/Users/etytel01/Documents/Acceleration/AlgoComparisons/Planar Experiment/encoder_8.dat'

    plt.ion()

    imu = IMU()

    imu.calibrate(calibfilename)
    imu.get_inertial_coords(calibfilename)
    imu.get_world_coordinates(calibfilename)

    imu.load(filename, resamplefreq=200.0)

    imu.filter(nsamp=10, method='running')

    with h5py.File(filename, 'r') as h5file:
        enc0 = np.array(h5file['/data/Encoder'])
        t0 = np.array(h5file['/data/t'])
    t0 /= 1000.0

    t = imu.t
    enc = interpolate.interp1d(t0, enc0)(t)

    # fig, ax = plt.subplots()
    # ax.plot(t, enc)

    orient_ekf1 = copy(imu.get_orientation(method='ekf', lCa=(0.0, -0.3, -0.7)))
    accd1 = copy(imu.accdyn)
    orient_ekf2 = copy(imu.get_orientation(method='ekf', lCa=(3.0, 3.0, 3.0)))
    accd2 = copy(imu.accdyn)
    orient_ekf3 = copy(imu.get_orientation(method='ekf', lCa=(0.0, 0.0, 0.0)))
    accd3 = copy(imu.accdyn)
    orient_ekf4 = copy(imu.get_orientation(method='ekf', lCa=(-5.0, -5.0, -5.0)))
    accd4 = copy(imu.accdyn)

    orient_mad = imu.get_orientation(method='madgwick')
    orient_gyro = imu.get_orientation(method='integrate_gyro')

    enc -= enc[0]

    fig, ax = plt.subplots()
    ax.plot(t, enc, label='encoder')
    ax.plot(imu.t, np.rad2deg(orient_ekf1[:, 0]), label='ekf')
    ax.plot(imu.t, np.rad2deg(orient_mad[:, 0]), label='madgwick')
    ax.plot(imu.t, np.rad2deg(orient_gyro[:, 0]), label='gyro')
    ax.legend()

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(t, enc, label='encoder')
    ax[0].plot(t, np.rad2deg(orient_ekf2[:,0]), label='Ca=1000')
    ax[0].plot(t, np.rad2deg(orient_ekf3[:,0]), label='Ca=1')
    ax[0].plot(t, np.rad2deg(orient_ekf4[:,0]), label='Ca=0.01')

    ax[1].plot(t, np.rad2deg(orient_ekf2[:,0])-enc, label='Ca=1000')
    ax[1].plot(t, np.rad2deg(orient_ekf3[:,0])-enc, label='Ca=1')
    ax[1].plot(t, np.rad2deg(orient_ekf4[:,0])-enc, label='Ca=0.01')
    ax[1].legend()

    fig, ax = plt.subplots()
    ax.plot(t, accd2[:, 0], label='Ca=1000')
    ax.plot(t, accd3[:, 0], label='Ca=1')
    ax.plot(t, accd4[:, 0], label='Ca=0.01')
    ax.legend()

    plt.show(block=True)

if __name__ == "__main__":
    main()
