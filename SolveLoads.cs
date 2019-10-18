using System;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace StructuralCalculator
{
    class PointLoad
    {
        public readonly double _p;
        public readonly double _L;
        public readonly double _a;
        public readonly double _b;
        public readonly double _rl;
        public readonly double _rr;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c4;

        public PointLoad(double p, double a, double L)
        {
            _p = p;
            _a = a;
            _L = L;
            _b = _L - _a;
            _rl = _p * _b / _L;
            _rr = _p * _a / _L;
            _c4 = (-1.0 * _rl * Math.Pow(_a, 3) / 3.0) - (_rr * Math.Pow(_a, 3) / 3.0) + (_rr * _L * Math.Pow(_a, 2) / 2.0);
            _c2 = -1.0 / _L * (_c4 + (_rr * Math.Pow(_L, 3) / 3.0));
            _c1 = (-1.0 * _rr * Math.Pow(_a, 2) / 2.0) - (_rl * Math.Pow(_a, 2) / 2.0) + (_rr * _L * _a) + _c2;
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    if (x[i] == 0 && _a == 0)
                    {
                        v[i] = 0.0;
                    }
                    else
                    {
                        v[i] = _rl;
                    }
                }
                else
                {
                    v[i] = -1.0 * _rr;
                }
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = _rl * x[i];
                }
                else
                {
                    m[i] = (-1.0 * _rr * x[i]) + (_rr * _L);
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = (_rl * Math.Pow(x[i], 2) / 2.0) + _c1;
                }
                else
                {
                    eis[i] = (-1.0 * _rr * Math.Pow(x[i], 2) / 2.0) + (_rr * _L * x[i]) + _c2;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = (_rl * Math.Pow(x[i], 3) / 6.0) + (_c1 * x[i]);
                }
                else
                {
                    eid[i] = (-1 * _rr * Math.Pow(x[i], 3) / 6.0) + (_rr * _L * Math.Pow(x[i], 2) / 2.0) + (_c2 * x[i]) + _c4;
                }
            }
            return eid;
        }

        public double Vx(double x)
        {
            double vx;
            if (x <= _a)
            {
                if (x == 0 && _a == 0)
                {
                    vx = 0.0;
                }
                else
                {
                    vx = _rl;
                }
            }
            else
            {
                vx = -1.0 * _rr;
            }
            return vx;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = _rl * x;
            }
            else
            {
                mx = (-1.0 * _rr * x) + (_rr * _L);
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = (_rl * Math.Pow(x, 2) / 2.0) + _c1;
            }
            else
            {
                eisx = (-1.0 * _rr * Math.Pow(x, 2) / 2.0) + (_rr * _L * x) + _c2;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = (_rl * Math.Pow(x, 3) / 6.0) + (_c1 * x);
            }
            else
            {
                eidx = (-1.0 * _rr * Math.Pow(x, 3) / 6.0) + (_rr * _L * Math.Pow(x, 2) / 2.0) + (_c2 * x) + _c4;
            }
            return eidx;
        }
    }

    class PointMoment
    {
        public readonly double _ma;
        public readonly double _L;
        public readonly double _a;
        public readonly double _rl;
        public readonly double _rr;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;

        public PointMoment(double ma, double a, double L)
        {
            _ma = ma;
            _a = a;
            _L = L;
            _rr = _ma / _L;
            _rl = -1.0 * _rr;
            _c2 = -1.0 / _L * ((_ma * Math.Pow(_a, 2)) - (0.5 * _ma * Math.Pow(_a, 2)) + (_rl * (Math.Pow(_L, 3) / 6.0)) + (0.5 * _ma * Math.Pow(_L, 2)));
            _c1 = _ma * _a + _c2;
            _c3 = 0.0;
            _c4 = (-1.0 * _rl * Math.Pow(_L, 3) / 6.0) - (0.5 * _ma * Math.Pow(_L, 2)) - (_c2 * _L);
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                v[i] = _rl;
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    if (x[i] == 0 && _a == 0)
                    {
                        m[i] = _ma;
                    }
                    else if (x[i] == _L && _a == _L)
                    {
                        m[i] = -1.0 * _ma;
                    }
                    else
                    {
                        m[i] = _rl * x[i];
                    }
                }
                else
                {
                    m[i] = (_rl * x[i]) + _ma;
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = (0.5 * _rl * Math.Pow(x[i], 2)) + _c1;
                }
                else
                {
                    eis[i] = (0.5 * _rl * Math.Pow(x[i], 2)) + (_ma * x[i]) + _c2;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = (1 / 6.0 * _rl * Math.Pow(x[i], 3)) + (_c1 * x[i]) + _c3;
                }
                else
                {
                    eid[i] = 1 / 6.0 * _rl * Math.Pow(x[i], 3) + (0.5 * _ma * Math.Pow(x[i], 2)) + (_c2 * x[i]) + _c4;
                }
            }
            return eid;
        }

        public double Vx()
        {
            return _rl;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                if (x == 0 && _a == 0)
                {
                    mx = _ma;
                }
                else if (x == _L && _a == _L)
                {
                    mx = -1.0 * _ma;
                }
                else
                {
                    mx = _rl * x;
                }
            }
            else
            {
                mx = (_rl * x) + _ma;
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = (0.5 * _rl * Math.Pow(x, 2)) + _c1;
            }
            else
            {
                eisx = (0.5 * _rl * Math.Pow(x, 2)) + (_ma * x) + _c2;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = (1 / 6.0 * _rl * Math.Pow(x, 3)) + (_c1 * x) + _c3;
            }
            else
            {
                eidx = 1 / 6.0 * _rl * Math.Pow(x, 3) + (0.5 * _ma * Math.Pow(x, 2)) + (_c2 * x) + _c4;
            }
            return eidx;
        }
    }

    class UniformDistributedLoad
    {
        public readonly double _w1;
        public readonly double _L;
        public readonly double _a;
        public readonly double _b;
        public readonly double _c;
        public readonly double _rl;
        public readonly double _rr;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;
        public readonly double _c5;
        public readonly double _c6;
        public readonly double _c7;
        public readonly double _c8;
        public readonly double _c9;

        public UniformDistributedLoad(double w1, double a, double b, double L)
        {
            _w1 = w1;
            _a = a;
            _L = L;
            _b = b;
            _c = _b - _a;
            _rl = (_w1 * _c) - (_w1 * _c * (_a + (_c / 2.0)) / _L);
            _rr = _w1 * _c * (_a + (_c / 2.0)) / _L;
            _c1 = 0.0;
            _c2 = -1.0 * _w1 * Math.Pow(_a, 2) / 2.0;
            _c3 = _rr * _L;
            _c7 = 0.0;
            _c8 = (-1.0 * _c1 * Math.Pow(_a, 2) / 2.0) + (_c2 * Math.Pow(_a, 2) / 2.0) + (5.0 * _w1 * Math.Pow(_a, 4) / 24.0) + _c7;
            _c9 = (-1.0 * _rl * Math.Pow(_b, 3) / 3.0) - (_rr * Math.Pow(_b, 3) / 3.0) + (_w1 * Math.Pow(_b, 4) / 8.0) - (_w1 * _a * Math.Pow(_b, 3) / 3.0) - (_c2 * Math.Pow(_b, 2) / 2.0) + (_c3 * Math.Pow(_b, 2) / 2.0) + _c8;
            _c6 = (_rr * Math.Pow(_L, 2) / 6.0) - (_c3 * _L / 2.0) - (_c9 / _L);
            _c5 = (-1.0 * _rl * Math.Pow(_b, 2) / 2.0) + (_w1 * Math.Pow(_b, 3) / 6.0) - (_w1 * _a * Math.Pow(_b, 2) / 2.0) - (_rr * Math.Pow(_b, 2) / 2.0) + (_c3 * _b) - (_c2 * _b) + _c6;
            _c4 = (_w1 * Math.Pow(_a, 3) / 3.0) + (_c2 * _a) + _c5 - (_c1 * _a);
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    v[i] = _rl;
                }
                else if (x[i] <= _b)
                {
                    v[i] = _rl - (_w1 * (x[i] - _a));
                }
                else
                {
                    v[i] = -1 * _rr;
                }
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = (_rl * x[i]) + _c1;
                }
                else if (x[i] <= _b)
                {
                    m[i] = (_rl * x[i]) - (_w1 * Math.Pow(x[i], 2) / 2) + (_w1 * _a * x[i]) + _c2;
                }
                else
                {
                    m[i] = (-1 * _rr * x[i]) + _c3;
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = (_rl * Math.Pow(x[i], 2) / 2.0) + (_c1 * x[i]) + _c4;
                }
                else if (x[i] <= _b)
                {
                    eis[i] = (_rl * Math.Pow(x[i], 2) / 2.0) - (_w1 * Math.Pow(x[i], 3) / 6.0) + (_w1 * _a * Math.Pow(x[i], 2) / 2.0) + (_c2 * x[i]) + _c5;
                }
                else
                {
                    eis[i] = (-1.0 * _rr * Math.Pow(x[i], 2) / 2.0) + (_c3 * x[i]) + _c6;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = (_rl * Math.Pow(x[i], 3) / 6.0) + (_c1 * Math.Pow(x[i], 2) / 2.0) + (_c4 * x[i]) + _c7;
                }
                else if (x[i] <= _b)
                {
                    eid[i] = (_rl * Math.Pow(x[i], 3) / 6.0) - (_w1 * Math.Pow(x[i], 4) / 24.0) + (_w1 * _a * Math.Pow(x[i], 3) / 6.0) + (_c2 * Math.Pow(x[i], 2) / 2.0) + (_c5 * x[i]) + _c8;
                }
                else
                {
                    eid[i] = (-1.0 * _rr * Math.Pow(x[i], 3) / 6.0) + (_c3 * Math.Pow(x[i], 2) / 2.0) + (_c6 * x[i]) + _c9;
                }
            }
            return eid;
        }

        public double Vx(double x)
        {
            double vx;
            if (x <= _a)
            {
                vx = _rl;
            }
            else if (x <= _b)
            {
                vx = _rl - (_w1 * (x - _a));
            }
            else
            {
                vx = -1.0 * _rr;
            }
            return vx;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = (_rl * x) + _c1;
            }
            else if (x <= _b)
            {
                mx = (_rl * x) - (_w1 * Math.Pow(x, 2) / 2.0) + (_w1 * _a * x) + _c2;
            }
            else
            {
                mx = (-1.0 * _rr * x) + _c3;
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = (_rl * Math.Pow(x, 2) / 2.0) + (_c1 * x) + _c4;
            }
            else if (x <= _b)
            {
                eisx = (_rl * Math.Pow(x, 2) / 2.0) - (_w1 * Math.Pow(x, 3) / 6.0) + (_w1 * _a * Math.Pow(x, 2) / 2.0) + (_c2 * x) + _c5;
            }
            else
            {
                eisx = (-1.0 * _rr * Math.Pow(x, 2) / 2.0) + (_c3 * x) + _c6;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = (_rl * Math.Pow(x, 3) / 6.0) + (_c1 * Math.Pow(x, 2) / 2.0) + (_c4 * x) + _c7;
            }
            else if (x <= _b)
            {
                eidx = (_rl * Math.Pow(x, 3) / 6.0) - (_w1 * Math.Pow(x, 4) / 24.0) + (_w1 * _a * Math.Pow(x, 3) / 6.0) + (_c2 * Math.Pow(x, 2) / 2.0) + (_c5 * x) + _c8;
            }
            else
            {
                eidx = (-1.0 * _rr * Math.Pow(x, 3) / 6.0) + (_c3 * Math.Pow(x, 2) / 2.0) + (_c6 * x) + _c9;
            }
            return eidx;
        }
    }

    class RightCantNoLoad
    {
        public readonly double _slope;
        public readonly double _rl;
        public readonly double _ml;

        public RightCantNoLoad(double slope)
        {
            _slope = slope;
            _rl = 0.0;
            _ml = 0.0;
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                eis[i] = _slope;
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                eid[i] = _slope * x[i];
            }
            return eid;
        }

        public double Vx()
        {
            return 0.0;
        }

        public double Mx()
        {
            return 0.0;
        }

        public double Eisx()
        {
            return _slope;
        }

        public double Eidx(double x)
        {
            return _slope * x;
        }
    }

    class RightCantPointLoad
    {
        public readonly double _p;
        public readonly double _L;
        public readonly double _Lb;
        public readonly double _a;
        public readonly double _b;
        public readonly double _rl;
        public readonly double _ml;
        public readonly PointMoment _backSpan;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;

        public RightCantPointLoad(double p, double a, double L, double Lb)
        {
            _p = p;
            _a = a;
            _L = L;
            _Lb = Lb;
            _b = _L - _a;
            _rl = _p;
            _ml = -1.0 * _p * _a;
            if (Lb == 0)
            {
                _c1 = 0.0;
            }
            else
            {
                _backSpan = new PointMoment(-1.0 * _ml, _Lb, _Lb);
                _c1 = _backSpan.Eisx(_Lb);
            }
            _c2 = 0.0;
            _c3 = 0.5 * _rl * Math.Pow(_a, 2) + _ml * _a + _c1;
            _c4 = -1.0 * _c3 * _a + 1.0 / 6.0 * _rl * Math.Pow(_a, 3) + 0.5 * _ml * Math.Pow(_a, 2) + _c1 * _a + _c2;
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    if (x[i] == 0 && _a == 0)
                    {
                        v[i] = 0.0;
                    }
                    else
                    {
                        v[i] = _p;
                    }
                }
                else
                {
                    v[i] = 0.0;
                }
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = _rl * x[i] + _ml;
                }
                else
                {
                    m[i] = 0.0;
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = 0.5 * _rl * Math.Pow(x[i], 2) + _ml * x[i] + _c1;
                }
                else
                {
                    eis[i] = _c3;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = 1.0 / 6.0 * _rl * Math.Pow(x[i], 3) + 0.5 * _ml * Math.Pow(x[i], 2) + _c1 * x[i] + _c2;
                }
                else
                {
                    eid[i] = _c3 * x[i] + _c4;
                }
            }
            return eid;
        }

        public double Vx(double x)
        {
            double vx;
            if (x <= _a)
            {
                if (x == 0 && _a == 0)
                {
                    vx = 0.0;
                }
                else
                {
                    vx = _p;
                }
            }
            else
            {
                vx = 0.0;
            }
            return vx;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = _rl * x + _ml;
            }
            else
            {
                mx = 0.0;
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = 0.5 * _rl * Math.Pow(x, 2) + _ml * x + _c1;
            }
            else
            {
                eisx = _c3;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = 1.0 / 6.0 * _rl * Math.Pow(x, 3) + 0.5 * _ml * Math.Pow(x, 2) + _c1 * x + _c2;
            }
            else
            {
                eidx = _c3 * x + _c4;
            }
            return eidx;
        }
    }

    class RightCantPointMoment
    {
        public readonly double _ma;
        public readonly double _L;
        public readonly double _Lb;
        public readonly double _a;
        public readonly double _b;
        public readonly double _rl;
        public readonly double _ml;
        public readonly PointMoment _backSpan;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;

        public RightCantPointMoment(double ma, double a, double L, double Lb)
        {
            _ma = ma;
            _a = a;
            _L = L;
            _Lb = Lb;
            _b = _L - _a;
            _rl = 0.0;
            _ml = -1.0 * _ma;
            if (Lb == 0)
            {
                _c1 = 0.0;
            }
            else
            {
                _backSpan = new PointMoment(-1.0 * _ml, _Lb, _Lb);
                _c1 = _backSpan.Eisx(_Lb);
            }
            _c2 = 0.0;
            _c3 = _ml * _a + _c1;
            _c4 = 0.5 * _ml * Math.Pow(_a, 2) + _c1 * _a + _c2 - _c3 * _a;
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                v[i] = 0.0;
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = _ml;
                }
                else
                {
                    m[i] = 0.0;
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = _ml * x[i] + _c1;
                }
                else
                {
                    eis[i] = _c3;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = 0.5 * _ml * Math.Pow(x[i], 2) + _c1 * x[i] + _c2;
                }
                else
                {
                    eid[i] = _c3 * x[i] + _c4;
                }
            }
            return eid;
        }

        public double Vx()
        {
            return 0.0;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = _ml;
            }
            else
            {
                mx = 0.0;
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = _ml * x + _c1;
            }
            else
            {
                eisx = _c3;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = 0.5 * _ml * Math.Pow(x, 2) + _c1 * x + _c2;
            }
            else
            {
                eidx = _c3 * x + _c4;
            }
            return eidx;
        }
    }

    class RightCantUniformDistributedLoad
    {
        public readonly double _w1;
        public readonly double _L;
        public readonly double _Lb;
        public readonly double _a;
        public readonly double _b;
        public readonly double _c;
        public readonly double _wTot;
        public readonly double _rl;
        public readonly double _ml;
        public readonly PointMoment _backSpan;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;
        public readonly double _c5;
        public readonly double _c6;

        public RightCantUniformDistributedLoad(double w1, double a, double b, double L, double Lb)
        {
            _w1 = w1;
            _L = L;
            _Lb = Lb;
            _a = a;
            _b = b;
            _c = _b - _a;
            _wTot = _w1 * _c;
            _rl = _wTot;
            _ml = -1.0 * _wTot * (_b - (_c / 2.0));
            if (Lb == 0)
            {
                _c1 = 0.0;
            }
            else
            {
                _backSpan = new PointMoment(-1.0 * _ml, _Lb, _Lb);
                _c1 = _backSpan.Eisx(_Lb);
            }
            _c2 = 0.0;
            _c3 = _c1;
            _c4 = _c1 * _a + _c2 - _c3 * a;
            _c5 = 0.5 * _wTot * Math.Pow(_b, 2) + _ml * _b - 1.0 / 6.0 * _w1 * Math.Pow(_b - _a, 3) + _c3;
            _c6 = 1.0 / 6.0 * _wTot * Math.Pow(_b, 3) + 0.5 * _ml * Math.Pow(_b, 2) - 1.0 / 24.0 * _w1 * Math.Pow(_b - _a, 4) + _c3 * _b + _c4 - _c5 * _b;
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    v[i] = _rl;
                }
                else if (x[i] <= _b)
                {
                    v[i] = _rl - _w1 * (x[i] - _a);
                }
                else
                {
                    v[i] = 0.0;
                }
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = _rl * x[i] + _ml;
                }
                else if (x[i] <= _b)
                {
                    m[i] = _rl * x[i] + _ml - (_w1 * (x[i] - _a) * ((x[i] - _a) / 2.0));
                }
                else
                {
                    m[i] = 0.0;
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = 0.5 * _rl * Math.Pow(x[i], 2) + _ml * x[i] + _c1;
                }
                else if (x[i] <= _b)
                {
                    eis[i] = 0.5 * _rl * Math.Pow(x[i], 2) + _ml * x[i] - (1.0 / 6.0 * _w1 * Math.Pow(x[i] - _a, 3)) + _c3;
                }
                else
                {
                    eis[i] = _c5;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = 1.0 / 6.0 * _rl * Math.Pow(x[i], 2) + 0.5 * _ml * Math.Pow(x[i], 2) + _c1 * x[i] + _c2;
                }
                else if (x[i] <= _b)
                {
                    eid[i] = 1.0 / 6.0 * _rl * Math.Pow(x[i], 3) + 0.5 * _ml * Math.Pow(x[i], 2) - 1.0 / 24.0 * _w1 * Math.Pow(x[i] - _a, 4) + _c3 * x[i] + _c4;
                }
                else
                {
                    eid[i] = _c5 * x[i] + _c6;
                }
            }
            return eid;
        }

        public double Vx(double x)
        {
            double vx;
            if (x <= _a)
            {
                vx = _wTot;
            }
            else if (x <= _b)
            {
                vx = _wTot - _w1 * (x - _a);
            }
            else
            {
                vx = 0.0;
            }
            return vx;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = _rl * x + _ml;
            }
            else if (x <= _b)
            {
                mx = _rl * x + _ml - (_w1 * (x - _a) * ((x - _a) / 2.0));
            }
            else
            {
                mx = 0.0;
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = 0.5 * _rl * Math.Pow(x, 2) + _ml * x + _c1;
            }
            else if (x <= _b)
            {
                eisx = 0.5 * _rl * Math.Pow(x, 2) + _ml * x - (1.0 / 6.0 * _w1 * Math.Pow(x - _a, 3)) + _c3;
            }
            else
            {
                eisx = _c5;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = 1.0 / 6.0 * _rl * Math.Pow(x, 2) + 0.5 * _ml * Math.Pow(x, 2) + _c1 * x + _c2;
            }
            else if (x <= _b)
            {
                eidx = 1.0 / 6.0 * _rl * Math.Pow(x, 3) + 0.5 * _ml * Math.Pow(x, 2) - 1.0 / 24.0 * _w1 * Math.Pow(x - _a, 4) + _c3 * x + _c4;
            }
            else
            {
                eidx = _c5 * x + _c6;
            }
            return eidx;
        }
    }

    class LeftCantNoLoad
    {
        public readonly double _L;
        public readonly double _slope;
        public readonly double _rr;
        public readonly double _mr;
        public readonly double _c1;
        public readonly double _c2;

        public LeftCantNoLoad(double slope, double L)
        {
            _slope = slope;
            _L = L;
            _rr = 0.0;
            _mr = 0.0;
            _c1 = _slope;
            _c2 = -1.0 * _c1 * _L;
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                eis[i] = _c1;
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                eid[i] = _c1 * x[i] + _c2;
            }
            return eid;
        }

        public double Vx()
        {
            return 0.0;
        }

        public double Mx()
        {
            return 0.0;
        }

        public double Eisx()
        {
            return _c1;
        }

        public double Eidx(double x)
        {
            return _c1 * x + _c2;
        }
    }

    class LeftCantPointLoad
    {
        public readonly double _p;
        public readonly double _L;
        public readonly double _Lb;
        public readonly double _a;
        public readonly double _b;
        public readonly double _rr;
        public readonly double _mr;
        public readonly PointMoment _backSpan;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;

        public LeftCantPointLoad(double p, double a, double L, double Lb)
        {
            _p = p;
            _a = a;
            _L = L;
            _Lb = Lb;
            _b = _L - _a;
            _rr = _p;
            _mr = -1.0 * _p * _b;
            if (Lb == 0)
            {
                _c3 = 0.0 + (0.5 * _p * Math.Pow(_L - _a, 2));
            }
            else
            {
                _backSpan = new PointMoment(_mr, 0.0, _Lb);
                _c3 = _backSpan.Eisx(0.0) + (0.5 * _p * Math.Pow(_L - _a, 2));
            }
            _c4 = (1 / 6.0 * _p * Math.Pow(_L - _a, 3)) - (_c3 * _L);
            _c1 = _c3;
            _c2 = (_c3 * _a) + _c4 - (_c1 * _a);
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    v[i] = 0.0;
                }
                else
                {
                    v[i] = -1.0 * _p;
                }
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = 0.0;
                }
                else
                {
                    m[i] = -1.0 * _p * (x[i] - _a);
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = _c1;
                }
                else
                {
                    eis[i] = (-0.5 * _p * Math.Pow(x[i] - _a, 2)) + _c3;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = _c1 * x[i] + _c2;
                }
                else
                {
                    eid[i] = -1.0 / 6.0 * _p * Math.Pow(x[i] - _a, 3) + _c3 * x[i] + _c4;
                }
            }
            return eid;
        }

        public double Vx(double x)
        {
            double vx;
            if (x <= _a)
            {
                vx = 0.0;
            }
            else
            {
                vx = -1.0 * _p;
            }
            return vx;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = 0.0;
            }
            else
            {
                mx = -1.0 * _p * (x - _a);
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = _c1;
            }
            else
            {
                eisx = (-0.5 * _p * Math.Pow(x - _a, 2)) + _c3;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = _c1 * x + _c2;
            }
            else
            {
                eidx = -1.0 / 6.0 * _p * Math.Pow(x - _a, 3) + _c3 * x + _c4;
            }
            return eidx;
        }
    }

    class LeftCantPointMoment
    {
        public readonly double _ma;
        public readonly double _L;
        public readonly double _Lb;
        public readonly double _a;
        public readonly double _rr;
        public readonly double _mr;
        public readonly PointMoment _backSpan;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;

        public LeftCantPointMoment(double ma, double a, double L, double Lb)
        {
            _ma = ma;
            _a = a;
            _L = L;
            _Lb = Lb;
            _rr = 0.0;
            _mr = _ma;
            if (Lb == 0)
            {
                _c3 = 0.0 - (_ma * _L);
            }
            else
            {
                _backSpan = new PointMoment(_mr, 0.0, _Lb);
                _c3 = _backSpan.Eisx(0.0) - (_ma * _L);
            }
            _c4 = (-0.5 * _ma * Math.Pow(_L, 2)) - _c3 * _L;
            _c1 = (1.0 * _ma * _a) + _c3;
            _c2 = 0.5 * _ma * Math.Pow(_a, 2) + _c3 * _a + _c4 - _c1 * _a;
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                v[i] = 0.0;
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = 0.0;
                }
                else
                {
                    m[i] = _ma;
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = _c1;
                }
                else
                {
                    eis[i] = (_ma * x[i]) + _c3;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = _c1 * x[i] + _c2;
                }
                else
                {
                    eid[i] = 0.5 * _ma * Math.Pow(x[i], 2) + _c3 * x[i] + _c4;
                }
            }
            return eid;
        }

        public double Vx()
        {
            return 0.0;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = 0.0;
            }
            else
            {
                mx = _ma;
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = _c1;
            }
            else
            {
                eisx = (_ma * x) + _c3;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = _c1 * x + _c2;
            }
            else
            {
                eidx = 0.5 * _ma * Math.Pow(x, 2) + _c3 * x + _c4;
            }
            return eidx;
        }
    }

    class LeftCantUniformDistributedLoad
    {
        public readonly double _w1;
        public readonly double _L;
        public readonly double _Lb;
        public readonly double _a;
        public readonly double _b;
        public readonly double _c;
        public readonly double _wTot;
        public readonly double _rr;
        public readonly double _mr;
        public readonly PointMoment _backSpan;
        public readonly double _c1;
        public readonly double _c2;
        public readonly double _c3;
        public readonly double _c4;
        public readonly double _c5;
        public readonly double _c6;

        public LeftCantUniformDistributedLoad(double w1, double a, double b, double L, double Lb)
        {
            _w1 = w1;
            _L = L;
            _Lb = Lb;
            _a = a;
            _b = b;
            _c = _b - _a;
            _wTot = _w1 * _c;
            _rr = _wTot;
            _mr = -1.0 * _wTot * (_L - (_a + (_c / 2.0)));
            if (Lb == 0)
            {
                _c5 = 0.0 + (0.5 * _wTot * Math.Pow(_L - (_a + (0.5 * _c)), 2));
            }
            else
            {
                _backSpan = new PointMoment(_mr, 0.0, _Lb);
                _c5 = _backSpan.Eisx(0.0) + (0.5 * _wTot * Math.Pow(_L - (_a + (0.5 * _c)), 2));
            }
            _c6 = (1.0 / 6.0 * _wTot * Math.Pow(_L - (_a + (0.5 * _c)), 3)) - (_c5 * _L);
            _c3 = ((-0.5) * _wTot * Math.Pow(_b - (_a + (0.5 * _c)), 2)) + _c5 + (1.0 / 6.0 * _w1 * Math.Pow(b - a, 3));
            _c1 = _c3;
            _c4 = (-1.0 / 6.0 * _wTot * Math.Pow(_b - (_a + (0.5 * _c)), 3)) + (_c5 * _b) + _c6 + (1.0 / 24.0 * _w1 * Math.Pow(_b - _a, 4)) - (_c3 * _b);
            _c2 = (_c3 * _a) + _c4 - (_c1 * _a);
        }

        public Vector<double> V(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> v = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    v[i] = 0.0;
                }
                else if (x[i] <= _b)
                {
                    v[i] = -1.0 * _w1 * (x[i] - _a);
                }
                else
                {
                    v[i] = -1.0 * _wTot;
                }
            }
            return v;
        }

        public Vector<double> M(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> m = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    m[i] = 0.0;
                }
                else if (x[i] <= _b)
                {
                    m[i] = -0.5 * _w1 * Math.Pow(x[i] - _a, 2);
                }
                else
                {
                    m[i] = -1.0 * _wTot * (x[i] - (_a + (0.5 * _c)));
                }
            }
            return m;
        }

        public Vector<double> Eis(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eis = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eis[i] = _c1;
                }
                else if (x[i] <= _b)
                {
                    eis[i] = -1.0 / 6.0 * _w1 * Math.Pow(x[i] - _a, 3) + _c3;
                }
                else
                {
                    eis[i] = (-0.5 * _wTot * Math.Pow(x[i] - (_a + (0.5 * _c)), 2)) + _c5;
                }
            }
            return eis;
        }

        public Vector<double> Eid(Vector<double> x)
        {
            int segments = x.Count;
            Vector<double> eid = Vector<double>.Build.Dense(segments);
            foreach (int i in Enumerable.Range(0, segments))
            {
                if (x[i] <= _a)
                {
                    eid[i] = _c1 * x[i] + _c2;
                }
                else if (x[i] <= _b)
                {
                    eid[i] = -1.0 / 24.0 * _w1 * Math.Pow(x[i] - _a, 4) + _c3 * x[i] + _c4;
                }
                else
                {
                    eid[i] = (-1.0 / 6.0 * _wTot * Math.Pow(x[i] - (_a + (0.5 * _c)), 3)) + _c5 * x[i] + _c6;
                }
            }
            return eid;
        }

        public double Vx(double x)
        {
            double vx;
            if (x <= _a)
            {
                vx = 0.0;
            }
            else if (x <= _b)
            {
                vx = -1.0 * _w1 * (x - _a);
            }
            else
            {
                vx = -1.0 * _wTot;
            }
            return vx;
        }

        public double Mx(double x)
        {
            double mx;
            if (x <= _a)
            {
                mx = 0.0;
            }
            else if (x <= _b)
            {
                mx = -0.5 * _w1 * Math.Pow(x - _a, 2);
            }
            else
            {
                mx = -1.0 * _wTot * (x - (_a + (0.5 * _c)));
            }
            return mx;
        }

        public double Eisx(double x)
        {
            double eisx;
            if (x <= _a)
            {
                eisx = _c1;
            }
            else if (x <= _b)
            {
                eisx = -1.0 / 6.0 * _w1 * Math.Pow(x - _a, 3) + _c3;
            }
            else
            {
                eisx = (-0.5 * _wTot * Math.Pow(x - (_a + (0.5 * _c)), 2)) + _c5;
            }
            return eisx;
        }

        public double Eidx(double x)
        {
            double eidx;
            if (x <= _a)
            {
                eidx = _c1 * x + _c2;
            }
            else if (x <= _b)
            {
                eidx = -1.0 / 24.0 * _w1 * Math.Pow(x - _a, 4) + _c3 * x + _c4;
            }
            else
            {
                eidx = (-1.0 / 6.0 * _wTot * Math.Pow(x - (_a + (0.5 * _c)), 3)) + _c5 * x + _c6;
            }
            return eidx;
        }
    }
}
