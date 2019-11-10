* Solving the Griffith problem using the XFEM *
* The XFEM C++ library developped by *
* S. Bordas C. Dunant N.V.Phu T.Q.Tri R. Duddu *
* ---------------------------------------------------------

ExportDisplacement 1 *
ExportStress 1 *
ExportStrain 1 *

Material 2 
1 class ElasticMaterial E 200.000000000000 n  0.300000000000 t 1.0 *
2 class ElasticMaterial E 100.000000000000 n  0.300000000000 t 1.0 *
**

TimeIntegrationScheme *
1 class Static *
**

Load 3 
1 class BoundaryCondition loadTimeFunction 1 conditions 1 d 0.00000000000000  * 
2 class NodalLoad   loadTimeFunction 1 components  2  0.00000000000000  0.60000000000000  * 
3 class NodalLoad   loadTimeFunction 1 components  2  0.00000000000000  0.30000000000000  * 
**

LoadTimeFunction 1 
1 class ConstantFunction f(t) 1.000000 *
**

Node 108 
1 nDofs 2 coord 2   0.000000000000   0.000000000000 bcOnDof1 1  bcOnDof2 1  * 
2 nDofs 2 coord 2   0.600000000000   0.000000000000 bcOnDof1 1  bcOnDof2 1  * 
3 nDofs 2 coord 2   1.200000000000   0.000000000000 bcOnDof1 1  bcOnDof2 1  * 
4 nDofs 2 coord 2   1.800000000000   0.000000000000 bcOnDof1 1  bcOnDof2 1  * 
5 nDofs 2 coord 2   2.400000000000   0.000000000000 bcOnDof1 1  bcOnDof2 1  * 
6 nDofs 2 coord 2   3.000000000000   0.000000000000 bcOnDof1 1  bcOnDof2 1  * 
7 nDofs 2 coord 2   0.000000000000   0.529411764706 * 
8 nDofs 2 coord 2   0.600000000000   0.529411764706 * 
9 nDofs 2 coord 2   1.200000000000   0.529411764706 * 
10 nDofs 2 coord 2   1.800000000000   0.529411764706 * 
11 nDofs 2 coord 2   2.400000000000   0.529411764706 * 
12 nDofs 2 coord 2   3.000000000000   0.529411764706 * 
13 nDofs 2 coord 2   0.000000000000   1.058823529412 * 
14 nDofs 2 coord 2   0.600000000000   1.058823529412 * 
15 nDofs 2 coord 2   1.200000000000   1.058823529412 * 
16 nDofs 2 coord 2   1.800000000000   1.058823529412 * 
17 nDofs 2 coord 2   2.400000000000   1.058823529412 * 
18 nDofs 2 coord 2   3.000000000000   1.058823529412 * 
19 nDofs 2 coord 2   0.000000000000   1.588235294118 * 
20 nDofs 2 coord 2   0.600000000000   1.588235294118 * 
21 nDofs 2 coord 2   1.200000000000   1.588235294118 * 
22 nDofs 2 coord 2   1.800000000000   1.588235294118 * 
23 nDofs 2 coord 2   2.400000000000   1.588235294118 * 
24 nDofs 2 coord 2   3.000000000000   1.588235294118 * 
25 nDofs 2 coord 2   0.000000000000   2.117647058824 * 
26 nDofs 2 coord 2   0.600000000000   2.117647058824 * 
27 nDofs 2 coord 2   1.200000000000   2.117647058824 * 
28 nDofs 2 coord 2   1.800000000000   2.117647058824 * 
29 nDofs 2 coord 2   2.400000000000   2.117647058824 * 
30 nDofs 2 coord 2   3.000000000000   2.117647058824 * 
31 nDofs 2 coord 2   0.000000000000   2.647058823529 * 
32 nDofs 2 coord 2   0.600000000000   2.647058823529 * 
33 nDofs 2 coord 2   1.200000000000   2.647058823529 * 
34 nDofs 2 coord 2   1.800000000000   2.647058823529 * 
35 nDofs 2 coord 2   2.400000000000   2.647058823529 * 
36 nDofs 2 coord 2   3.000000000000   2.647058823529 * 
37 nDofs 2 coord 2   0.000000000000   3.176470588235 * 
38 nDofs 2 coord 2   0.600000000000   3.176470588235 * 
39 nDofs 2 coord 2   1.200000000000   3.176470588235 * 
40 nDofs 2 coord 2   1.800000000000   3.176470588235 * 
41 nDofs 2 coord 2   2.400000000000   3.176470588235 * 
42 nDofs 2 coord 2   3.000000000000   3.176470588235 * 
43 nDofs 2 coord 2   0.000000000000   3.705882352941 * 
44 nDofs 2 coord 2   0.600000000000   3.705882352941 * 
45 nDofs 2 coord 2   1.200000000000   3.705882352941 * 
46 nDofs 2 coord 2   1.800000000000   3.705882352941 * 
47 nDofs 2 coord 2   2.400000000000   3.705882352941 * 
48 nDofs 2 coord 2   3.000000000000   3.705882352941 * 
49 nDofs 2 coord 2   0.000000000000   4.235294117647 * 
50 nDofs 2 coord 2   0.600000000000   4.235294117647 * 
51 nDofs 2 coord 2   1.200000000000   4.235294117647 * 
52 nDofs 2 coord 2   1.800000000000   4.235294117647 * 
53 nDofs 2 coord 2   2.400000000000   4.235294117647 * 
54 nDofs 2 coord 2   3.000000000000   4.235294117647 * 
55 nDofs 2 coord 2   0.000000000000   4.764705882353 * 
56 nDofs 2 coord 2   0.600000000000   4.764705882353 * 
57 nDofs 2 coord 2   1.200000000000   4.764705882353 * 
58 nDofs 2 coord 2   1.800000000000   4.764705882353 * 
59 nDofs 2 coord 2   2.400000000000   4.764705882353 * 
60 nDofs 2 coord 2   3.000000000000   4.764705882353 * 
61 nDofs 2 coord 2   0.000000000000   5.294117647059 * 
62 nDofs 2 coord 2   0.600000000000   5.294117647059 * 
63 nDofs 2 coord 2   1.200000000000   5.294117647059 * 
64 nDofs 2 coord 2   1.800000000000   5.294117647059 * 
65 nDofs 2 coord 2   2.400000000000   5.294117647059 * 
66 nDofs 2 coord 2   3.000000000000   5.294117647059 * 
67 nDofs 2 coord 2   0.000000000000   5.823529411765 * 
68 nDofs 2 coord 2   0.600000000000   5.823529411765 * 
69 nDofs 2 coord 2   1.200000000000   5.823529411765 * 
70 nDofs 2 coord 2   1.800000000000   5.823529411765 * 
71 nDofs 2 coord 2   2.400000000000   5.823529411765 * 
72 nDofs 2 coord 2   3.000000000000   5.823529411765 * 
73 nDofs 2 coord 2   0.000000000000   6.352941176471 * 
74 nDofs 2 coord 2   0.600000000000   6.352941176471 * 
75 nDofs 2 coord 2   1.200000000000   6.352941176471 * 
76 nDofs 2 coord 2   1.800000000000   6.352941176471 * 
77 nDofs 2 coord 2   2.400000000000   6.352941176471 * 
78 nDofs 2 coord 2   3.000000000000   6.352941176471 * 
79 nDofs 2 coord 2   0.000000000000   6.882352941176 * 
80 nDofs 2 coord 2   0.600000000000   6.882352941176 * 
81 nDofs 2 coord 2   1.200000000000   6.882352941176 * 
82 nDofs 2 coord 2   1.800000000000   6.882352941176 * 
83 nDofs 2 coord 2   2.400000000000   6.882352941176 * 
84 nDofs 2 coord 2   3.000000000000   6.882352941176 * 
85 nDofs 2 coord 2   0.000000000000   7.411764705882 * 
86 nDofs 2 coord 2   0.600000000000   7.411764705882 * 
87 nDofs 2 coord 2   1.200000000000   7.411764705882 * 
88 nDofs 2 coord 2   1.800000000000   7.411764705882 * 
89 nDofs 2 coord 2   2.400000000000   7.411764705882 * 
90 nDofs 2 coord 2   3.000000000000   7.411764705882 * 
91 nDofs 2 coord 2   0.000000000000   7.941176470588 * 
92 nDofs 2 coord 2   0.600000000000   7.941176470588 * 
93 nDofs 2 coord 2   1.200000000000   7.941176470588 * 
94 nDofs 2 coord 2   1.800000000000   7.941176470588 * 
95 nDofs 2 coord 2   2.400000000000   7.941176470588 * 
96 nDofs 2 coord 2   3.000000000000   7.941176470588 * 
97 nDofs 2 coord 2   0.000000000000   8.470588235294 * 
98 nDofs 2 coord 2   0.600000000000   8.470588235294 * 
99 nDofs 2 coord 2   1.200000000000   8.470588235294 * 
100 nDofs 2 coord 2   1.800000000000   8.470588235294 * 
101 nDofs 2 coord 2   2.400000000000   8.470588235294 * 
102 nDofs 2 coord 2   3.000000000000   8.470588235294 * 
103 nDofs 2 coord 2   0.000000000000   9.000000000000 loads 1 2  * 
104 nDofs 2 coord 2   0.600000000000   9.000000000000 loads 1 2  * 
105 nDofs 2 coord 2   1.200000000000   9.000000000000 loads 1 2  * 
106 nDofs 2 coord 2   1.800000000000   9.000000000000 loads 1 2  * 
107 nDofs 2 coord 2   2.400000000000   9.000000000000 loads 1 2  * 
108 nDofs 2 coord 2   3.000000000000   9.000000000000 loads 1 2  * 
**

Element 85 
1 class Q4U nodes 1 2 8 7  mat 1  * 
2 class Q4U nodes 2 3 9 8  mat 1  * 
3 class Q4U nodes 3 4 10 9  mat 1  * 
4 class Q4U nodes 4 5 11 10  mat 1  * 
5 class Q4U nodes 5 6 12 11  mat 1  * 
6 class Q4U nodes 7 8 14 13  mat 1  * 
7 class Q4U nodes 8 9 15 14  mat 1  * 
8 class Q4U nodes 9 10 16 15  mat 1  * 
9 class Q4U nodes 10 11 17 16  mat 1  * 
10 class Q4U nodes 11 12 18 17  mat 1  * 
11 class Q4U nodes 13 14 20 19  mat 1  * 
12 class Q4U nodes 14 15 21 20  mat 1  * 
13 class Q4U nodes 15 16 22 21  mat 1  * 
14 class Q4U nodes 16 17 23 22  mat 1  * 
15 class Q4U nodes 17 18 24 23  mat 1  * 
16 class Q4U nodes 19 20 26 25  mat 1  * 
17 class Q4U nodes 20 21 27 26  mat 1  * 
18 class Q4U nodes 21 22 28 27  mat 1  * 
19 class Q4U nodes 22 23 29 28  mat 1  * 
20 class Q4U nodes 23 24 30 29  mat 1  * 
21 class Q4U nodes 25 26 32 31  mat 1  * 
22 class Q4U nodes 26 27 33 32  mat 1  * 
23 class Q4U nodes 27 28 34 33  mat 1  * 
24 class Q4U nodes 28 29 35 34  mat 1  * 
25 class Q4U nodes 29 30 36 35  mat 1  * 
26 class Q4U nodes 31 32 38 37  mat 1  * 
27 class Q4U nodes 32 33 39 38  mat 1  * 
28 class Q4U nodes 33 34 40 39  mat 1  * 
29 class Q4U nodes 34 35 41 40  mat 1  * 
30 class Q4U nodes 35 36 42 41  mat 1  * 
31 class Q4U nodes 37 38 44 43  mat 1  * 
32 class Q4U nodes 38 39 45 44  mat 1  * 
33 class Q4U nodes 39 40 46 45  mat 1  * 
34 class Q4U nodes 40 41 47 46  mat 1  * 
35 class Q4U nodes 41 42 48 47  mat 1  * 
36 class Q4U nodes 43 44 50 49  mat 1  * 
37 class Q4U nodes 44 45 51 50  mat 1  * 
38 class Q4U nodes 45 46 52 51  mat 1  * 
39 class Q4U nodes 46 47 53 52  mat 1  * 
40 class Q4U nodes 47 48 54 53  mat 1  * 
41 class Q4U nodes 49 50 56 55  mat 1  * 
42 class Q4U nodes 50 51 57 56  mat 1  * 
43 class Q4U nodes 51 52 58 57  mat 1  * 
44 class Q4U nodes 52 53 59 58  mat 1  * 
45 class Q4U nodes 53 54 60 59  mat 1  * 
46 class Q4U nodes 55 56 62 61  mat 1  * 
47 class Q4U nodes 56 57 63 62  mat 1  * 
48 class Q4U nodes 57 58 64 63  mat 1  * 
49 class Q4U nodes 58 59 65 64  mat 1  * 
50 class Q4U nodes 59 60 66 65  mat 1  * 
51 class Q4U nodes 61 62 68 67  mat 1  * 
52 class Q4U nodes 62 63 69 68  mat 1  * 
53 class Q4U nodes 63 64 70 69  mat 1  * 
54 class Q4U nodes 64 65 71 70  mat 1  * 
55 class Q4U nodes 65 66 72 71  mat 1  * 
56 class Q4U nodes 67 68 74 73  mat 1  * 
57 class Q4U nodes 68 69 75 74  mat 1  * 
58 class Q4U nodes 69 70 76 75  mat 1  * 
59 class Q4U nodes 70 71 77 76  mat 1  * 
60 class Q4U nodes 71 72 78 77  mat 1  * 
61 class Q4U nodes 73 74 80 79  mat 1  * 
62 class Q4U nodes 74 75 81 80  mat 1  * 
63 class Q4U nodes 75 76 82 81  mat 1  * 
64 class Q4U nodes 76 77 83 82  mat 1  * 
65 class Q4U nodes 77 78 84 83  mat 1  * 
66 class Q4U nodes 79 80 86 85  mat 1  * 
67 class Q4U nodes 80 81 87 86  mat 1  * 
68 class Q4U nodes 81 82 88 87  mat 1  * 
69 class Q4U nodes 82 83 89 88  mat 1  * 
70 class Q4U nodes 83 84 90 89  mat 1  * 
71 class Q4U nodes 85 86 92 91  mat 1  * 
72 class Q4U nodes 86 87 93 92  mat 1  * 
73 class Q4U nodes 87 88 94 93  mat 1  * 
74 class Q4U nodes 88 89 95 94  mat 1  * 
75 class Q4U nodes 89 90 96 95  mat 1  * 
76 class Q4U nodes 91 92 98 97  mat 1  * 
77 class Q4U nodes 92 93 99 98  mat 1  * 
78 class Q4U nodes 93 94 100 99  mat 1  * 
79 class Q4U nodes 94 95 101 100  mat 1  * 
80 class Q4U nodes 95 96 102 101  mat 1  * 
81 class Q4U nodes 97 98 104 103  mat 1  * 
82 class Q4U nodes 98 99 105 104  mat 1  * 
83 class Q4U nodes 99 100 106 105  mat 1  * 
84 class Q4U nodes 100 101 107 106  mat 1  * 
85 class Q4U nodes 101 102 108 107  mat 1  * 
**

EnrichmentItem 3 
1 class CrackInterior  geometry 1 myTips 1 2 EnrichmentFuncs 1  1  enrichScheme 3 * 
2 class CrackTip  Type BiMatElast Mat 2  1 2 geometry 3 EnrichmentFuncs  12  2 3 4 5 6 7 8 9 10 11 12 13 enrichScheme 1 domainIntFac 5.3f 2.500000e+00 class CrackTip  Type BiMatElast Mat 2  1 2 geometry  * 
3 class MaterialInterface  geometry 4 Mat 2 1 2 EnrichmentFuncs 1 14  enrichScheme 3 * 
**

GeometryEntity 6 
1 class PiecewiseLinear  numOfVertices 2 vertices 2 3 geoDescription 1 *
2 class Vertex  coord 2 -0.010000000  4.500000000 geoDescription 1 * 
3 class Vertex  coord 2  1.500000000  4.500000000 geoDescription 1 * 
4 class PiecewiseLinear  numOfVertices 2 vertices 5 6 geoDescription 1 * 
5 class Vertex  coord 2  -0.01 4.50 geoDescription 1 * 
6 class Vertex  coord 2  3.00 4.50 geoDescription 1 * 
**

EnrichmentFunction 14 
1 class DiscontinuousField  * 
2 class BiMatCrackAsymp1 * 
3 class BiMatCrackAsymp2 * 
4 class BiMatCrackAsymp3 * 
5 class BiMatCrackAsymp4 * 
6 class BiMatCrackAsymp5 * 
7 class BiMatCrackAsymp6 * 
8 class BiMatCrackAsymp7 * 
9 class BiMatCrackAsymp8 * 
10 class BiMatCrackAsymp9 * 
11 class BiMatCrackAsymp10 * 
12 class BiMatCrackAsymp11 * 
13 class BiMatCrackAsymp12 * 
14 class AbsSignedDistance  * 
**

NLSolver *
1 class NewtonRaphson n 100 t 1e-5 c 1 *
**

TimeStep 1 
1 dt 1.0 *
**

