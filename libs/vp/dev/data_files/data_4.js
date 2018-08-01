var pts = 'M26.4028,29.895 L24.2703,54.2699 L21.1965,89.4031 L43.7381,91.3752 L120.553,98.0956 L122.416,76.7961 L125.759,38.5875 L93.9578,35.8053 L93.1577,35.7353 L85.9114,35.1014 L83.7556,59.7427 L77.6073,59.2048 L69.0207,58.4535 L64.0949,58.0226 L65.2072,45.3088 L65.4052,43.0457 L65.5974,40.8491 L65.9235,37.1215 L66.2508,33.3813 Z'; var polys = {
	"0" : "M26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 L65.5974,40.8491 L65.4052,43.0457 L65.2072,45.3088 L64.0949,58.0226 L117.428,97.8223 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 Z",
	"1" : "M24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 L65.5974,40.8491 L65.4052,43.0457 L65.2072,45.3088 L64.0949,58.0226 L123.568,63.6268 L122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 L21.1965,89.4031 Z",
	"10" : "M83.7556,59.7427 L85.9114,35.1014 L93.1577,35.7353 L93.9578,35.8053 L125.759,38.5875 L122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 L21.1965,89.4031 L24.247,54.5363 Z",
	"11" : "M77.6073,59.2048 L123.603,63.2289 L122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 L21.1965,89.4031 L24.247,54.5363 Z",
	"12" : "M69.0207,58.4535 L123.603,63.2289 L122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 L21.1965,89.4031 L24.247,54.5363 Z",
	"13" : "M64.0949,58.0226 L123.603,63.2289 L122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 L65.5974,40.8491 L65.4052,43.0457 L65.2072,45.3088 Z",
	"14" : "M65.2072,45.3088 L61.0445,92.8893 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 L65.5974,40.8491 L65.4052,43.0457 Z",
	"15" : "M65.4052,43.0457 L61.0445,92.8893 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 L65.5974,40.8491 Z",
	"16" : "M65.5974,40.8491 L61.0445,92.8893 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 Z",
	"17" : "M65.9235,37.1215 L61.0445,92.8893 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 Z",
	"18" : "M66.2508,33.3813 L61.0445,92.8893 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 Z",
	"2" : "M21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 L65.5974,40.8491 L65.4052,43.0457 L65.2072,45.3088 L64.0949,58.0226 L69.0207,58.4535 L77.6073,59.2048 L83.7556,59.7427 L125.646,39.8818 L122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 Z",
	"3" : "M43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.4028,29.895 L66.2508,33.3813 L65.9235,37.1215 L65.5974,40.8491 L65.4052,43.0457 L65.2072,45.3088 L64.0949,58.0226 L69.0207,58.4535 L77.6073,59.2048 L83.7556,59.7427 L112.037,37.387 L125.759,38.5875 L122.416,76.7961 L120.553,98.0956 Z",
	"4" : "M120.553,98.0956 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L26.2896,31.1888 L64.0949,58.0226 L69.0207,58.4535 L77.6073,59.2048 L83.7556,59.7427 L85.9114,35.1014 L93.1577,35.7353 L93.9578,35.8053 L125.759,38.5875 L122.416,76.7961 Z",
	"5" : "M122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 L21.1965,89.4031 L24.2703,54.2699 L25.0418,45.4514 L64.0949,58.0226 L69.0207,58.4535 L77.6073,59.2048 L83.7556,59.7427 L85.9114,35.1014 L93.1577,35.7353 L93.9578,35.8053 L125.759,38.5875 Z",
	"6" : "M125.759,38.5875 L122.416,76.7961 L120.553,98.0956 L43.7381,91.3752 L24.3222,89.6765 L83.7556,59.7427 L85.9114,35.1014 L93.1577,35.7353 L93.9578,35.8053 Z",
	"7" : "M93.9578,35.8053 L125.759,38.5875 L122.416,76.7961 L120.553,98.0956 L69.3196,93.6133 L83.7556,59.7427 L85.9114,35.1014 L93.1577,35.7353 Z",
	"8" : "M93.1577,35.7353 L93.9578,35.8053 L125.759,38.5876 L122.416,76.7961 L120.553,98.0956 L70.4518,93.7124 L83.7556,59.7427 L85.9114,35.1014 Z",
	"9" : "M85.9114,35.1014 L93.1577,35.7353 L93.9578,35.8053 L125.759,38.5876 L122.416,76.7961 L120.553,98.0956 L80.7051,94.6094 Z"
};