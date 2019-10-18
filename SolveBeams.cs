using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace StructuralCalculator
{
    /// <summary>
    /// 
    /// Implementation of the Theory of Three Moments, https://en.wikipedia.org/wiki/Theorem_of_three_moments
    /// 
    ///     Inputs:
    /// 
    ///         beam_depth = beam depth as a double -- Expected Units: in
    /// 
    ///         beam_width = beam width as a double -- Expected Units: in
    /// 
    ///         beam_modulus = modulus of elasticity assumed to be constant over all spans as a double -- Expected Units: psi
    /// 
    ///         beam_momentofinteria = moment of inertia assumed to be constant over all spans as a double -- Expected Units: in^4
    /// 
    ///         beam_spans = span lengths as a List<double> -- Expected Units: in -- Example: { 120.00, 120.00 }
    /// 
    ///         beam_loads_raw = applied loads as a List<List<object>> -- Examples:
    /// 
    ///             Point Loads:
    ///             { P, P, a, a, 'POINT', span }
    ///                 P = Load -- Expected Units: lbs
    ///                 a = Load location from left support -- Expected Units: in
    ///                 span = integer indicating which span the load is in -- 0 for first span
    ///                       P+
    ///                 ___a__|_____
    ///                 ^           ^
    ///                 
    ///             Point Moments:
    ///             { M, M, a, a, 'POINT_MOMENT', span }
    ///                 M = Moment -- Expected Units: in-lbs
    ///                 a = Moment location from left support -- Expected Units: in
    ///                 span = integer indicating which span the load is in -- 0 for first span
    ///                       --->M+
    ///                 ___a__|___|___
    ///                 ^     {---     ^
    ///                 
    ///             Uniform loads:
    ///             { w, w, a, b, 'UDL', span }
    ///                 w = Load -- Expected Units: lbs per in
    ///                 a = Load start location from left support -- Expected Units: in
    ///                 b = Load end location from left support -- Expected Units: in
    ///                 span = integer indicating which span the load is in -- 0 for first span
    ///                      ___w+__
    ///                 ___a_|_____|_____
    ///                 ^          |     ^
    ///                 |----b-----|
    ///                 
    ///         cantilever = cantilever presence as a boolean -- first span in beam_spans is analyzed as cantilever if True
    ///         
    ///         segments = number of segments to create for span in beam_spans as an integer
    /// 
    /// </summary>

    class BeamSolver
    {
        private double _beamDepth;
        private double _beamWidth;
        private double _beamModulus;
        private double _beamMomentOfInertia;
        private List<double> _beamSpans;
        private List<List<object>> _beamLoadsRaw;
        private bool _cantilever;
        private int _segments;
        private int _numberOfSpans;
        private List<double> _cumulativeSumOfSpans;
        private Matrix<double> _segmentsNonCumulative = null;
        private Matrix<double> _segmentsCumulative = null;
        private double? _cantileverMoments = null;
        private Matrix<double> _loadMoments = null;
        private Matrix<double> _deflectionDemand = null;
        private Matrix<double> _bendingMomentDemand = null;
        private Matrix<double> _shearDemand = null;
        private Matrix<double> _slopeOrRotation = null;
        private Vector<double> _supportMoments = null;
        private Vector<double> _supportReactions = null;

        public BeamSolver(double beamDepth, double beamWidth, double beamModulus, double beamMomentOfInertia, List<double> beamSpans, List<List<object>> beamLoadsRaw, bool cantilever, int segments)
        {
            _beamDepth = beamDepth;
            _beamWidth = beamWidth;
            _beamModulus = beamModulus;
            _beamMomentOfInertia = beamMomentOfInertia;
            _beamSpans = beamSpans;
            _beamLoadsRaw = beamLoadsRaw;
            _cantilever = cantilever;
            _segments = segments;
            _numberOfSpans = _beamSpans.Count;
            double sum = 0.0;
            _cumulativeSumOfSpans = _beamSpans
                .Select(w => sum += w).ToList();
        }

        public Matrix<double> BeamSegments
        {
            get
            {
                if (_segmentsCumulative is null && _segmentsNonCumulative is null)
                {
                    GetSegmentsCumulative();
                    GetSegmentsNonCumulative();
                }
                return _segmentsCumulative;
            }
        }

        public List<Matrix<double>> Demand
        {
            get
            {
                if (_deflectionDemand is null && _bendingMomentDemand is null && _shearDemand is null && _slopeOrRotation is null)
                {
                    if (_segmentsCumulative is null && _segmentsNonCumulative is null)
                    {
                        GetSegmentsCumulative();
                        GetSegmentsNonCumulative();
                    }
                    if (_cantileverMoments is null)
                    {
                        GetCantileverMoments();
                    }
                    if (_loadMoments is null)
                    {
                        GetLoadMoments();
                    }
                    if (_supportMoments is null)
                    {
                        GetSupportMoments();
                    }
                    GetDemand();
                }
                return new List<Matrix<double>> { _shearDemand, _bendingMomentDemand, _slopeOrRotation, _deflectionDemand };
            }
        }

        public Vector<double> SupportMoments
        {
            get
            {
                if (_supportMoments is null)
                {
                    if (_segmentsCumulative is null && _segmentsNonCumulative is null)
                    {
                        GetSegmentsCumulative();
                        GetSegmentsNonCumulative();
                    }
                    if (_cantileverMoments is null)
                    {
                        GetCantileverMoments();
                    }
                    if (_loadMoments is null)
                    {
                        GetLoadMoments();
                    }
                    GetSupportMoments();
                }
                return _supportMoments;
            }
        }

        public Vector<double> SupportReactions
        {
            get
            {
                if (_supportReactions is null)
                {
                    if (_segmentsCumulative is null && _segmentsNonCumulative is null)
                    {
                        GetSegmentsCumulative();
                        GetSegmentsNonCumulative();
                    }
                    if (_supportMoments is null)
                    {
                        if (_cantileverMoments is null)
                        {
                            GetCantileverMoments();
                        }
                        if (_loadMoments is null)
                        {
                            GetLoadMoments();
                        }
                        GetSupportMoments();
                    }
                    GetSupportReactions();
                }
                return _supportReactions;
            }
        }

        private void GetSegmentsNonCumulative()
        {
            Matrix<double> beamSegments = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                beamSegments[0, j] = 0.0;
                beamSegments[_segments, j] = _beamSpans[j];
                foreach (int i in Enumerable.Range(1, _segments - 1))
                {
                    beamSegments[i, j] = beamSegments[i - 1, j] + _beamSpans[j] / _segments;
                }
            }
            _segmentsNonCumulative = beamSegments; // used in the calculations
        }

        private void GetSegmentsCumulative()
        {
            Matrix<double> beamSegments = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                beamSegments[0, j] = 0.0;
                beamSegments[_segments, j] = _beamSpans[j];
                foreach (int i in Enumerable.Range(1, _segments - 1))
                {
                    beamSegments[i, j] = beamSegments[i - 1, j] + _beamSpans[j] / _segments;
                }
            }
            // Converts lengths to global rather than local ie span 2 x[0] now = span 1 x[-1] or length in lieu of 0
            foreach (int j in Enumerable.Range(1, _numberOfSpans - 1))
            {
                foreach (int i in Enumerable.Range(0, _segments + 1))
                {
                    beamSegments[i, j] += _cumulativeSumOfSpans[j - 1];
                }
            }
            beamSegments = beamSegments.Divide(12); // convert from inches to feet
            _segmentsCumulative = beamSegments;
        }

        private void GetCantileverMoments()
        {
            double cantileverMoment = 0.0;
            if (_cantilever is true)
            {
                foreach (List<object> loads in _beamLoadsRaw)
                {
                    if ((int)loads[5] == 0)
                    {
                        if ((string)loads[4] == "POINT")
                        {
                            LeftCantPointLoad load = new LeftCantPointLoad((double)loads[0], (double)loads[2], _beamSpans[0], 0);
                            cantileverMoment += load._mr;
                        }
                        else if ((string)loads[4] == "POINT_MOMENT")
                        {
                            LeftCantPointMoment load = new LeftCantPointMoment((double)loads[0], (double)loads[2], _beamSpans[0], 0);
                            cantileverMoment += load._mr;
                        }
                        else if ((string)loads[4] == "UDL")
                        {
                            LeftCantUniformDistributedLoad load = new LeftCantUniformDistributedLoad((double)loads[0], (double)loads[2], (double)loads[3], _beamSpans[0], 0);
                            cantileverMoment += load._mr;
                        }
                    }
                }
                _cantileverMoments = cantileverMoment;
            }
        }

        private void GetLoadMoments()
        {
            Matrix<double> segmentedLoadMoments = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            // Span as simple support
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                foreach (List<object> loads in _beamLoadsRaw)
                {
                    if ((int)loads[5] == j)
                    {
                        if ((string)loads[4] == "POINT")
                        {
                            PointLoad load = new PointLoad((double)loads[0], (double)loads[2], _beamSpans[j]);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(j));
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedLoadMoments[i, j] += moment[i];
                            }
                        }
                        else if ((string)loads[4] == "POINT_MOMENT")
                        {
                            PointMoment load = new PointMoment((double)loads[0], (double)loads[2], _beamSpans[j]);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(j));
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedLoadMoments[i, j] += moment[i];
                            }
                        }
                        else if ((string)loads[4] == "UDL")
                        {
                            UniformDistributedLoad load = new UniformDistributedLoad((double)loads[0], (double)loads[2], (double)loads[3], _beamSpans[j]);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(j));
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedLoadMoments[i, j] += moment[i];
                            }
                        }
                    }
                }
            }
            _loadMoments = segmentedLoadMoments;
        }

        private void GetDemand()
        {
            Matrix<double> segmentedDeflection = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            Matrix<double> segmentedCombinedMoments = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            Matrix<double> segmentedShear = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            Matrix<double> segmentedSlopeOrRotation = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            // Span as simple support Moment, Shears, and Reactions
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                foreach (List<object> loads in _beamLoadsRaw)
                {
                    if ((int)loads[5] == j)
                    {
                        if ((string)loads[4] == "POINT")
                        {
                            PointLoad load = new PointLoad((double)loads[0], (double)loads[2], _beamSpans[j]);
                            Vector<double> deflection = load.Eid(_segmentsNonCumulative.Column(j)).Divide(_beamMomentOfInertia * _beamModulus);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(j));
                            Vector<double> shear = load.V(_segmentsNonCumulative.Column(j));
                            Vector<double> slopeOrRotation = load.Eis(_segmentsNonCumulative.Column(j)).Divide(_beamMomentOfInertia * _beamModulus);
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedDeflection[i, j] += deflection[i];
                                segmentedCombinedMoments[i, j] += moment[i];
                                segmentedShear[i, j] += shear[i];
                                segmentedSlopeOrRotation[i, j] += slopeOrRotation[i];
                            }
                        }
                        else if ((string)loads[4] == "POINT_MOMENT")
                        {
                            PointMoment load = new PointMoment((double)loads[0], (double)loads[2], _beamSpans[j]);
                            Vector<double> deflection = load.Eid(_segmentsNonCumulative.Column(j)).Divide(_beamMomentOfInertia * _beamModulus);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(j));
                            Vector<double> shear = load.V(_segmentsNonCumulative.Column(j));
                            Vector<double> slopeOrRotation = load.Eis(_segmentsNonCumulative.Column(j)).Divide(_beamMomentOfInertia * _beamModulus);
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedDeflection[i, j] += deflection[i];
                                segmentedCombinedMoments[i, j] += moment[i];
                                segmentedShear[i, j] += shear[i];
                                segmentedSlopeOrRotation[i, j] += slopeOrRotation[i];
                            }
                        }
                        else if ((string)loads[4] == "UDL")
                        {
                            UniformDistributedLoad load = new UniformDistributedLoad((double)loads[0], (double)loads[2], (double)loads[3], _beamSpans[j]);
                            Vector<double> deflection = load.Eid(_segmentsNonCumulative.Column(j)).Divide(_beamMomentOfInertia * _beamModulus);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(j));
                            Vector<double> shear = load.V(_segmentsNonCumulative.Column(j));
                            Vector<double> slopeOrRotation = load.Eis(_segmentsNonCumulative.Column(j)).Divide(_beamMomentOfInertia * _beamModulus);
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedDeflection[i, j] += deflection[i];
                                segmentedCombinedMoments[i, j] += moment[i];
                                segmentedShear[i, j] += shear[i];
                                segmentedSlopeOrRotation[i, j] += slopeOrRotation[i];
                            }
                        }
                    }
                }
            }
            // Cantilever Moments, Shears, and Reactions
            if (_cantilever is true)
            {
                foreach (int i in Enumerable.Range(0, _segments + 1))
                {
                    segmentedDeflection[i, 0] = 0.0;
                    segmentedCombinedMoments[i, 0] = 0.0;
                    segmentedShear[i, 0] = 0.0;
                    segmentedSlopeOrRotation[i, 0] = 0.0;
                }
                foreach (List<object> loads in _beamLoadsRaw)
                {
                    if ((int)loads[5] == 0)
                    {
                        if ((string)loads[4] == "POINT")
                        {
                            LeftCantPointLoad load = new LeftCantPointLoad((double)loads[0], (double)loads[2], _beamSpans[0], 0);
                            Vector<double> deflection = load.Eid(_segmentsNonCumulative.Column(0)).Divide(_beamMomentOfInertia * _beamModulus);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(0));
                            Vector<double> shear = load.V(_segmentsNonCumulative.Column(0));
                            Vector<double> slopeOrRotation = load.Eis(_segmentsNonCumulative.Column(0)).Divide(_beamMomentOfInertia * _beamModulus);
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedDeflection[i, 0] += deflection[i];
                                segmentedCombinedMoments[i, 0] += moment[i];
                                segmentedShear[i, 0] += shear[i];
                                segmentedSlopeOrRotation[i, 0] += slopeOrRotation[i];
                            }
                        }
                        else if ((string)loads[4] == "POINT_MOMENT")
                        {
                            LeftCantPointMoment load = new LeftCantPointMoment((double)loads[0], (double)loads[2], _beamSpans[0], 0);
                            Vector<double> deflection = load.Eid(_segmentsNonCumulative.Column(0)).Divide(_beamMomentOfInertia * _beamModulus);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(0));
                            Vector<double> shear = load.V(_segmentsNonCumulative.Column(0));
                            Vector<double> slopeOrRotation = load.Eis(_segmentsNonCumulative.Column(0)).Divide(_beamMomentOfInertia * _beamModulus);
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedDeflection[i, 0] += deflection[i];
                                segmentedCombinedMoments[i, 0] += moment[i];
                                segmentedShear[i, 0] += shear[i];
                                segmentedSlopeOrRotation[i, 0] += slopeOrRotation[i];
                            }
                        }
                        else if ((string)loads[4] == "UDL")
                        {
                            LeftCantUniformDistributedLoad load = new LeftCantUniformDistributedLoad((double)loads[0], (double)loads[2], (double)loads[3], _beamSpans[0], 0);
                            Vector<double> deflection = load.Eid(_segmentsNonCumulative.Column(0)).Divide(_beamMomentOfInertia * _beamModulus);
                            Vector<double> moment = load.M(_segmentsNonCumulative.Column(0));
                            Vector<double> shear = load.V(_segmentsNonCumulative.Column(0));
                            Vector<double> slopeOrRotation = load.Eis(_segmentsNonCumulative.Column(0)).Divide(_beamMomentOfInertia * _beamModulus);
                            foreach (int i in Enumerable.Range(0, _segments + 1))
                            {
                                segmentedDeflection[i, 0] += deflection[i];
                                segmentedCombinedMoments[i, 0] += moment[i];
                                segmentedShear[i, 0] += shear[i];
                                segmentedSlopeOrRotation[i, 0] += slopeOrRotation[i];
                            }
                        }
                    }
                }
            }
            // Support Moment response on spans
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                if (j == 0 && _cantilever is true)
                {
                    continue;
                }
                else
                {
                    PointMoment load1 = new PointMoment(_supportMoments[j], 0, _beamSpans[j]);
                    PointMoment load2 = new PointMoment(-1.0 * _supportMoments[j + 1], _beamSpans[j], _beamSpans[j]);
                    Vector<double> deflection = load1.Eid(_segmentsNonCumulative.Column(j)).Add(load2.Eid(_segmentsNonCumulative.Column(j))).Divide(_beamMomentOfInertia * _beamModulus);
                    Vector<double> moment = load1.M(_segmentsNonCumulative.Column(j)).Add(load2.M(_segmentsNonCumulative.Column(j)));
                    Vector<double> shear = load1.V(_segmentsNonCumulative.Column(j)).Add(load2.V(_segmentsNonCumulative.Column(j)));
                    Vector<double> slopeOrRotation = load1.Eis(_segmentsNonCumulative.Column(j)).Add(load2.Eis(_segmentsNonCumulative.Column(j))).Divide(_beamMomentOfInertia * _beamModulus);
                    foreach (int i in Enumerable.Range(0, _segments + 1))
                    {
                        segmentedDeflection[i, j] += deflection[i];
                        segmentedCombinedMoments[i, j] += moment[i];
                        segmentedShear[i, j] += shear[i];
                        segmentedSlopeOrRotation[i, j] += slopeOrRotation[i];
                    }
                }
            }
            // Fixes cantilever slope and deflection for interior span end slopes
            if (_cantilever is true)
            {
                LeftCantNoLoad cantileverFix = new LeftCantNoLoad(segmentedSlopeOrRotation[0, 1], _beamSpans[0]);
                Vector<double> deflection = cantileverFix.Eid(_segmentsNonCumulative.Column(0));
                Vector<double> slopeOrRotation = cantileverFix.Eis(_segmentsNonCumulative.Column(0));
                foreach (int i in Enumerable.Range(0, _segments + 1))
                {
                    segmentedDeflection[i, 0] += deflection[i];
                    segmentedSlopeOrRotation[i, 0] += slopeOrRotation[i];
                }
            }
            // Convert units
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                foreach (int i in Enumerable.Range(0, _segments + 1))
                {
                    segmentedCombinedMoments[i, j] = segmentedCombinedMoments[i, j] / (12.0 * 1000.0); // in-lbs to ft-kips
                    segmentedShear[i, j] = segmentedShear[i, j] / 1000.0; // lbs to kips
                }
            }
            _deflectionDemand = segmentedDeflection;
            _bendingMomentDemand = segmentedCombinedMoments;
            _shearDemand = segmentedShear;
            _slopeOrRotation = segmentedSlopeOrRotation;
        }

        private void GetSupportMoments()
        {
            // Horizontal center of moment region
            Matrix<double> aXlXr = Matrix<double>.Build.Dense(3, _numberOfSpans);
            Matrix<double> mXx = Matrix<double>.Build.Dense(_segments + 1, _numberOfSpans);
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                foreach (int i in Enumerable.Range(0, _segments + 1))
                {
                    mXx[i, j] = _loadMoments[i, j] * _segmentsNonCumulative[i, j];
                }
                double A = Simps(_loadMoments.Column(j), _segmentsNonCumulative.Column(j));
                aXlXr[0, j] = A;
                if (A == 0)
                {
                    aXlXr[1, j] = 0.0;
                    aXlXr[2, j] = 0.0;
                }
                else
                {
                    double xL = 1.0 / A * Simps(mXx.Column(j), _segmentsNonCumulative.Column(j));
                    aXlXr[1, j] = xL;
                    aXlXr[2, j] = _beamSpans[j] - xL;
                }
            }
            Matrix<double> F = Matrix<double>.Build.Dense(_numberOfSpans + 1, _numberOfSpans + 1);
            foreach (int j in Enumerable.Range(1, _numberOfSpans - 1))
            {
                F[j, j - 1] = _beamSpans[j - 1] / _beamMomentOfInertia;
                F[j, j] = 2.0 * ((_beamSpans[j - 1] / _beamMomentOfInertia) + (_beamSpans[j] / _beamMomentOfInertia));
                F[j, j + 1] = _beamSpans[j] / _beamMomentOfInertia;
            }
            if (_cantilever is true)
            {
                F[0, 0] = 1.0;
                foreach (int j in Enumerable.Range(0, _numberOfSpans))
                {
                    F[1, j] = 0.0;
                }
                F[1, 1] = 1.0;
                F[_numberOfSpans, _numberOfSpans] = 1.0;
            }
            else
            {
                F[0, 0] = 1.0;
                F[_numberOfSpans, _numberOfSpans] = 1.0;
            }
            // Support settlement delta
            Vector<double> delta = Vector<double>.Build.Dense(_numberOfSpans + 1);
            foreach (int j in Enumerable.Range(1, _numberOfSpans - 1))
            {
                delta[j] = -6.0 * aXlXr[1, j - 1] * aXlXr[0, j - 1] / (_beamSpans[j - 1] * _beamMomentOfInertia) + -6.0 * aXlXr[2, j] * aXlXr[0, j] / (_beamSpans[j] * _beamMomentOfInertia);
            }
            if (_cantilever is true)
            {
                delta[0] = 0.0;
                delta[1] = (double)_cantileverMoments;
            }
            // Support moments
            _supportMoments = F.Inverse().Multiply(delta);
        }

        private void GetSupportReactions()
        {
            // Reactions at interior supports
            Vector<double> supportReactions = Vector<double>.Build.Dense(_numberOfSpans + 1);
            Matrix<double> rSpan = Matrix<double>.Build.Dense(2, _numberOfSpans);
            // Span as simple support
            foreach (int j in Enumerable.Range(0, _numberOfSpans))
            {
                foreach (List<object> loads in _beamLoadsRaw)
                {
                    if ((int)loads[5] == j)
                    {
                        if ((string)loads[4] == "POINT")
                        {
                            PointLoad load = new PointLoad((double)loads[0], (double)loads[2], _beamSpans[j]);
                            rSpan[0, j] += load._rl;
                            rSpan[1, j] += load._rr;
                        }
                        else if ((string)loads[4] == "POINT_MOMENT")
                        {
                            PointMoment load = new PointMoment((double)loads[0], (double)loads[2], _beamSpans[j]);
                            rSpan[0, j] += load._rl;
                            rSpan[1, j] += load._rr;

                        }
                        else if ((string)loads[4] == "UDL")
                        {
                            UniformDistributedLoad load = new UniformDistributedLoad((double)loads[0], (double)loads[2], (double)loads[3], _beamSpans[j]);
                            rSpan[0, j] += load._rl;
                            rSpan[1, j] += load._rr;
                        }
                    }
                }
            }
            foreach (int j in Enumerable.Range(0, _numberOfSpans + 1))
            {
                if (j == 0)
                {
                    supportReactions[j] = rSpan[0, j] - (_supportMoments[j] / _beamSpans[j]) + (_supportMoments[j + 1] / _beamSpans[j]);
                }
                else if (j == _numberOfSpans)
                {
                    supportReactions[j] = rSpan[1, j - 1] - (_supportMoments[j] / _beamSpans[j - 1]) + (_supportMoments[j - 1] / _beamSpans[j - 1]);
                }
                else if (j > 0 && j < _numberOfSpans)
                {
                    supportReactions[j] = rSpan[0, j] + rSpan[1, j - 1] - (_supportMoments[j] / _beamSpans[j]) - (_supportMoments[j] / _beamSpans[j - 1]) + (_supportMoments[j + 1] / _beamSpans[j]) + (_supportMoments[j - 1] / _beamSpans[j - 1]);
                }
            }
            _supportReactions = supportReactions;
        }

        private double BasicSimps(double[] y, int start, int stop, int dx, double[] x = null)
        {
            double result;
            int step = 2;
            double[] slice0 = y.Skip(start).Where((value, index) => (index % step == 0) && (index < stop)).ToArray();
            double[] slice1 = y.Skip(start + 1).Where((value, index) => index % step == 0 && (index < stop + 1)).ToArray();
            double[] slice2 = y.Skip(start + 2).Where((value, index) => index % step == 0 && (index < stop + 2)).ToArray();
            if (x is null)
            {
                result = slice0.Zip(slice1.Select(r => r * 4.0).ToArray(), (aX, aY) => aX + aY).ToArray()
                    .Zip(slice2, (aX, aY) => aX + aY).Select(r => r * dx / 3.0).ToArray()
                    .Sum();
            }
            else
            {
                double[] h = x.Skip(1).ToArray().Zip(x.Where((value, index) => index < x.Length).ToArray(), (aX, aY) => aX - aY).ToArray();
                double[] h0 = h.Skip(start).Where((value, index) => (index % step == 0) && (index < stop)).ToArray();
                double[] h1 = h.Skip(start + 1).Where((value, index) => index % step == 0 && (index < stop + 1)).ToArray();
                double[] hsum = h0.Zip(h1, (aX, aY) => aX + aY).ToArray();
                double[] hprod = h0.Zip(h1, (aX, aY) => aX * aY).ToArray();
                double[] hquot = h0.Zip(h1, (aX, aY) => aX / aY).ToArray();

                result = slice0.Zip(hquot.Select(r => 2.0 - (1.0 / r)).ToArray(), (aX, aY) => aX * aY).ToArray()
                    .Zip(slice1.Zip(hsum, (aX, aY) => aX * aY).ToArray().Zip(hsum, (aX, aY) => aX * aY).ToArray().Zip(hprod, (aX, aY) => aX / aY).ToArray(), (aX, aY) => aX + aY).ToArray()
                    .Zip(slice2.Zip(hquot.Select(r => 2.0 - r).ToArray(), (aX, aY) => aX * aY).ToArray(), (aX, aY) => aX + aY).ToArray()
                    .Zip(hsum.Select(r => r / 6.0).ToArray(), (aX, aY) => aX * aY).ToArray()
                    .Sum();
            }
            return result;
        }

        private double Simps(Vector<double> y, Vector<double> x = null, int dx = 1, string even = "avg")
        {
            double[] xA = x?.ToArray();
            double[] yA = y.ToArray();
            double firstDx = dx;
            double lastDx = dx;
            int N = yA.Length;
            int nd = yA.Rank;
            double result = 0.0;
            double val = 0.0;
            if (N % 2 == 0)
            {
                if (!new string[] { "avg", "last", "first" }.Contains(even))
                {
                    throw new ArgumentException("Parameter 'even' must be 'avg', 'last', or 'first'.", "even");
                }
                if (new string[] { "avg", "first" }.Contains(even))
                {
                    if (x != null)
                    {
                        lastDx = xA[N - 1] - xA[N - 2];
                    }
                    val += 0.5 * lastDx * (yA[N - 1] + yA[N - 2]);
                    result = BasicSimps(yA, 0, N - 3, dx, xA);
                }
                if (new string[] { "avg", "last" }.Contains(even))
                {
                    if (x != null)
                    {
                        firstDx = xA[1] - xA[0];
                    }
                    val += 0.5 * firstDx * (yA[1] + yA[0]);
                    result += BasicSimps(yA, 1, N - 2, dx, xA);
                }
                if (even == "avg")
                {
                    val /= 2.0;
                    result /= 2.0;
                }
                result += val;
            }
            else
            {
                result = BasicSimps(yA, 0, N - 2, dx, xA);
            }
            return result;
        }
    }
}
