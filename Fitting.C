#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

/*
 * Code fitting entered results from experiments
 * 												*/

Double_t fitFunction(Double_t *x, Double_t *par) {
  return par[0] + 0*par[1]*x[0] + 1*par[2]*x[0]*x[0];
  //return par[0] + 1*par[1]*x[0] + 0*par[2]*x[0]*x[0]+par[2]*x[0]*x[0]*x[0];
  //return log(x[0])/(par[0]*sqrt(x[0]*log(x[0])));  //Landau
  //return TMath::Binomial(par[0],x[0])*TMath::Power(par[1],x[0])*TMath::Power((1-par[1]),(par[0]-x[0])); //n=par[0], p=par[1], k=x[0]
  //return TMath::Binomiall(par[1], par[0], x[0]);
}

void FittingDemo() {
	
   const TString titlex = "Products of random Artin generators";
   const TString titley = "Occurence";

   TCanvas *c1 = new TCanvas("c1","Fitting",10,10,700,500);
   c1->SetFillColor(0);
   c1->SetFrameFillColor(0);
   c1->SetGrid();
 
   const int nBins = 17;
   
   //only success runtime N as datasets
   Double_t data1[nBins] 	= {0.234, 0.496, 0.90333, 1.425, 2.025, 2.59, 3.66, 4.42, 5.74, 7.4975, 8.25, 10.2525, 0, 13.68167, 0, 19.19, 23.51};
   Double_t data2[nBins]	= {0.29111, 0.5725, 0.925, 1.55167, 2.13167, 3.21857, 4.12286, 5.40857, 5.71, 8.664, 10.3525, 11.89556, 14.64333, 15.04667, 17.13, 19.678, 0}; //letzter Wert als Ausreißer weggelassen
   Double_t data3[nBins]	= {0.31556,  0.608,  1.06,  1.66625,  2.2225,  3.00833,  4.26714,  5.85625,  7.052,  9.045,  9.91667,  11.53167,  12.25,  17.24667,  18.42667,  21.44143,  23.33333};
   Double_t data4[nBins]	= {0.33375, 0.65875, 1.37143, 1.78333, 2.44714, 3.55167, 5.174, 5.77, 7.31833, 9.57, 11.48, 14.62, 14.07333, 15.572, 22.815, 23.32667, 24.6575};
   Double_t data5[nBins]	= {0.3375, 0.69333, 1.308, 1.93111, 3.58, 3.60857, 5.43833, 6.1075, 8.22, 10.49142, 10.37, 16.4775, 16.26, 22.21, 20.17, 25.82333, 31.89};
   Double_t data6[nBins]	= {0.405, 0.90714, 1.49889, 2.26125, 3.48, 4.918, 5.96286, 7.935, 10.15375, 12.468, 9.74333, 16.045, 23.03333, 24.97, 19.88, 27.57667, 35.63};
   Double_t data7[nBins]	= {0.46142, 0.94625, 1.82714, 2.73125, 3.13, 4.72429, 7.7375, 9.675, 8.68, 14.77, 6.75, 19.89667, 19.52667, 25.11, 28.98, 32.9625, 0};
   Double_t data8[nBins]	= {0.43778, 1.02556, 1.69286, 2.62833, 4.22, 4.6925, 8.08, 9.156, 10.875, 15.23, 11.61333, 21.07, 23.185, 27.74, 29.48, 38.74333, 0};
   Double_t data9[nBins]	= {0.49, 1.1, 1.74125, 2.815, 3.80286, 5.355, 6.944, 12.62, 12.8675, 16.546, 22.34, 19.195, 25.6, 33.99, 43.25, 42.36, 31.95};
   Double_t data10[nBins]	= {0.47666, 1.13571, 2.01125, 2.73, 4.0125, 6.068, 7.076, 10.25, 14.78333, 15.01667, 12.37, 21.895, 22.87667, 25.19, 0, 47.49, 51.81};
   Double_t data11[nBins]	= {0.4725, 1.115, 1.96143, 2.725, 3.8575, 6.396, 7.76333, 9.78, 13.82, 17.17, 19.82333, 0, 0, 0, 0, 40.8, 44.81};
   Double_t data12[nBins]	= {0.535, 1.01333, 1.9475, 2.97, 5.005, 5.65, 7.68, 11.73, 12.5125, 15.395, 19.125, 20.13, 33.55, 0, 40.4, 0, 62.08};
   Double_t data13[nBins]	= {0.545, 1.18, 1.9675, 2.90167, 4.2625, 6.33, 8.52667, 11.79, 13.17, 22.16, 16.925, 25.49, 26.405, 0, 34.63, 34.3, 43.59};
   Double_t data14[nBins]	= {0.54778, 1.18, 2.395, 3.96167, 5.68, 7.5875, 8.58, 12.68, 0, 20.06, 28.05, 0, 35.53, 38.7, 52.45, 55.4, 0};
   Double_t data15[nBins]	= {0.5725, 1.40571, 2.01667, 3.676, 5.5375, 9.405, 9.08, 12.7625, 11.3, 25.04, 25.21333, 35.96, 0, 0, 57.95, 44.87, 72.47};
   Double_t data16[nBins]	= {0.57, 1.34857, 2.416, 4.338, 6.33, 6.91, 11.3075, 11.54, 17.255, 16.85, 23.08, 24.59, 34.02, 39.84, 0, 0, 0};
   Double_t data17[nBins]	= {0.57714, 1.27667, 2.22333, 3.51286, 7.32333, 7.30333, 8.49, 13.565, 17.31, 0, 26.59, 30.4, 0, 32.87, 0, 47.45, 0};
   
   //runtime N as datasets
   //Double_t data1[nBins] 	= {0.234, 0.427, 0.903, 1.018, 1.511, 2.041, 2.903, 2.83, 3.675, 5.586, 5.715, 6.773, 5.161, 11.026, 7.272, 13.217, 9.882};
   //Double_t data2[nBins]	= {0.296, 0.542, 0.925, 1.428, 2.052, 2.766, 3.753, 4.223, 3.281, 5.998, 6.355, 11.547, 9.677, 11.138, 13.38, 13.284, 9.649}; 
   //Double_t data3[nBins]	= {0.312, 0.608, 0.814, 1.593, 1.955, 2.425, 3.757, 5.075, 5.357, 3.584, 7.946, 8.875, 8.7, 11.839, 10.745, 17.26, 18.589};
   //Double_t data4[nBins]	= {0.338, 0.608, 1.106, 1.65, 1.924, 2.713, 3.459, 4.97, 6.571, 7.652, 7.63, 12.03, 9.984, 11.31, 11.617, 17.408, 17.921};
   //Double_t data5[nBins]	= {0.335, 0.643, 1.308, 1.94, 2.474, 2.881, 4.618, 5.263, 6.257, 9.093, 8.136, 11.264, 9.715, 10.215, 11.685, 18.396, 13.935};
   //Double_t data6[nBins]	= {0.422, 0.869, 1.409, 2.325, 3.177, 4.426, 5.629, 5.448, 8.651, 10.128, 8.82, 10.619, 17.97, 14.784, 8.013, 18.834, 18.093};
   //Double_t data7[nBins]	= {0.459, 0.88, 1.777, 2.584, 2.991, 4.149, 5.702, 7.867, 8.122, 9.429, 9.029, 13.741, 15.908, 16.488, 13.898, 25.043, 26.181};
   //Double_t data8[nBins]	= {0.452, 0.986, 1.565, 2.608, 3.378, 3.896, 5.413, 7.458, 6.665, 12.283, 11.599, 13.42, 16.057, 18.878, 18.168, 22.317, 21.463};
   //Double_t data9[nBins]	= {0.49, 1.028, 1.843, 2.658, 3.731, 5.405, 6.565, 8.148, 8.651, 14.782, 9.649, 14.812, 19.997, 19.996, 17.709, 32.136, 20.448};
   //Double_t data10[nBins]	= {0.492, 1.093, 1.906, 2.547, 3.554, 4.332, 5.663, 6.673, 8.042, 13.845, 14.408, 15.894, 14.78, 17.836, 17.042, 25.517, 35.385};
   //Double_t data11[nBins]	= {0.445, 1.108, 1.961, 2.778, 3.663, 5.029, 7.335, 5.771, 8.424, 10.132, 15.339, 12.248, 14.476, 22.655, 27.814, 28.707, 26.541};
   //Double_t data12[nBins]	= {0.5, 0.981, 1.786, 2.235, 4.835, 4.564, 6.922, 8.924, 10.012, 12.889, 14.74, 15.312, 19.459, 17.691, 22.971, 31.028, 30.561};
   //Double_t data13[nBins]	= {0.544, 1.102, 1.776, 3.064, 3.211, 5.681, 5.768, 7.461, 9.954, 13.731, 13.312, 15.63, 20.077, 16.82, 20.21, 32.582, 27.509};
   //Double_t data14[nBins]	= {0.516, 1.396, 2.708, 3.365, 4.549, 5.89, 10.083, 10.335, 8.91, 10.483, 14.103, 18.757, 27.723, 28.894, 41.264, 25.119, 24.573};
   //Double_t data15[nBins]	= {0.578, 1.443, 2.194, 3.69, 5.216, 6.177, 6.815, 11.598, 14.186, 15.12, 19.455, 16.674, 19.783, 21.877, 34.129, 37.317, 47.155};
   //Double_t data16[nBins]	= {0.582, 1.336, 2.281, 3.773, 6.038, 5.363, 8.834, 8.392, 11.914, 13.835, 14.561, 31.981, 17.257, 27.014, 39.508, 22.57, 24.309};
   //Double_t data17[nBins]	= {0.534, 1.402, 2.16, 3.883, 5.46, 6.571, 8.487, 10.021, 11.999, 14.542, 18.374, 18.951, 24.403, 23.065, 28.893, 33.572, 40.226};
   
   //runtime |x| as datasets
   //Double_t data1[nBins] 	= {0.234, 0.296, 0.312, 0.338, 0.335, 0.422, 0.459, 0.452, 0.49, 0.492, 0.445, 0.5, 0.544, 0.516, 0.578, 0.582, 0.534};
   //Double_t data2[nBins]	= {0.427, 0.542, 0.608, 0.608, 0.643, 0.869, 0.88, 0.986, 1.028, 1.093, 1.108, 0.981, 1.102, 1.396, 1.443, 1.336, 1.402}; 
   //Double_t data3[nBins]	= {0.903, 0.925, 0.814, 1.106, 1.308, 1.409, 1.777, 1.565, 1.843, 1.906, 1.961, 1.786, 1.776, 2.708, 2.194, 2.281, 2.16};
   //Double_t data4[nBins]	= {1.018, 1.428, 1.593, 1.65, 1.94, 2.325, 2.584, 2.608, 2.658, 2.547, 2.778, 2.235, 3.064, 3.365, 3.69, 3.773, 3.883};
   //Double_t data5[nBins]	= {1.511, 2.052, 1.955, 1.924, 2.474, 3.177, 2.991, 3.378, 3.731, 3.554, 3.663, 4.835, 3.211, 4.549, 5.216, 6.038, 5.46};
   //Double_t data6[nBins]	= {2.041, 2.766, 2.425, 2.713, 2.881, 4.426, 4.149, 3.896, 5.405, 4.332, 5.029, 4.564, 5.681, 5.89, 6.177, 5.363, 6.571};
   //Double_t data7[nBins]	= {2.903, 3.753, 3.757, 3.459, 4.618, 5.629, 5.702, 5.413, 6.565, 5.663, 7.335, 6.922, 5.768, 10.083, 6.815, 8.834, 8.487};
   //Double_t data8[nBins]	= {2.83, 4.223, 5.075, 4.97, 5.263, 5.448, 7.867, 7.458, 8.148, 6.673, 5.771, 8.924, 7.461, 10.335, 11.598, 8.392, 10.021};
   //Double_t data9[nBins]	= {3.675, 3.281, 5.357, 6.571, 6.257, 8.651, 8.122, 6.665, 8.651, 8.042, 8.424, 10.012, 9.954, 8.91, 14.186, 11.914, 11.999};
   //Double_t data10[nBins]	= {5.586, 5.998, 3.584, 7.652, 9.093, 10.128, 9.429, 12.283, 14.782, 13.845, 10.132, 12.889, 13.731, 10.483, 15.12, 13.835, 14.542};
   //Double_t data11[nBins]	= {5.715, 6.355, 7.946, 7.63, 8.136, 8.82, 9.029, 11.599, 9.649, 14.408, 15.339, 14.74, 13.312, 14.103, 19.455, 14.561, 18.374};
   //Double_t data12[nBins]	= {6.773, 11.547, 8.875, 12.03, 11.264, 10.619, 13.741, 13.42, 14.812, 15.894, 12.248, 15.312, 15.63, 18.757, 16.674, 31.981, 18.951};
   //Double_t data13[nBins]	= {5.161, 9.677, 8.7, 9.984, 9.715, 17.97, 15.908, 16.057, 19.997, 14.78, 14.476, 19.459, 20.077, 27.723, 19.783, 17.257, 24.403};
   //Double_t data14[nBins]	= {11.026, 11.138, 11.839, 11.31, 10.215, 14.784, 16.488, 18.878, 19.996, 17.836, 22.655, 17.691, 16.82, 28.894, 21.877, 27.014, 23.065};
   //Double_t data15[nBins]	= {7.272, 13.38, 10.745, 11.617, 11.685, 8.013, 13.898, 18.168, 17.709, 17.042, 27.814, 22.971, 20.21, 41.264, 34.129, 39.508, 28.893};
   //Double_t data16[nBins]	= {13.217, 13.284, 17.26, 17.408, 18.396, 18.834, 25.043, 22.317, 32.136, 25.517, 28.707, 31.028, 32.582, 25.119, 37.317, 22.57, 33.572};
   //Double_t data17[nBins]	= {9.882, 9.649, 18.589, 17.921, 13.935, 18.093, 26.181, 21.463, 20.448, 35.385, 26.541, 30.561, 27.509, 24.573, 47.155, 24.309, 40.226};
   
   //probability of success N as datasets
   //Double_t data1[nBins] 	= {0.888, 0.798, 0.734, 0.654, 0.648, 0.558, 0.484, 0.46, 0.458, 0.414, 0.388, 0.324, 0.296, 0.28, 0.24, 0.222, 0.246};
   //Double_t data2[nBins]	= {0.896, 0.818, 0.78, 0.734, 0.712, 0.638, 0.622, 0.536, 0.536, 0.538, 0.48, 0.43, 0.45, 0.348, 0.33, 0.37, 0.298}; 
   //Double_t data3[nBins]	= {0.884, 0.85, 0.782, 0.726, 0.742, 0.63, 0.586, 0.58, 0.516, 0.502, 0.464, 0.456, 0.39, 0.368, 0.366, 0.376, 0.358};
   //Double_t data4[nBins]	= {0.882, 0.848, 0.736, 0.736, 0.696, 0.616, 0.616, 0.558, 0.546, 0.494, 0.434, 0.422, 0.41, 0.396, 0.332, 0.32, 0.334};
   //Double_t data5[nBins]	= {0.888, 0.83, 0.77, 0.7, 0.648, 0.598, 0.526, 0.498, 0.478, 0.466, 0.446, 0.404, 0.344, 0.336, 0.308, 0.244, 0.256};
   //Double_t data6[nBins]	= {0.862, 0.768, 0.702, 0.716, 0.622, 0.598, 0.524, 0.512, 0.466, 0.446, 0.376, 0.352, 0.332, 0.324, 0.31, 0.246, 0.246};
   //Double_t data7[nBins]	= {0.898, 0.802, 0.754, 0.676, 0.586, 0.528, 0.484, 0.458, 0.424, 0.382, 0.328, 0.304, 0.296, 0.25, 0.244, 0.234, 0.216};
   //Double_t data8[nBins]	= {0.856, 0.76, 0.71, 0.634, 0.566, 0.514, 0.504, 0.406, 0.376, 0.34, 0.322, 0.272, 0.2, 0.23, 0.2, 0.184, 0.168};
   //Double_t data9[nBins]	= {0.864, 0.78, 0.678, 0.632, 0.552, 0.488, 0.406, 0.392, 0.35, 0.308, 0.264, 0.228, 0.24, 0.196, 0.18, 0.164, 0.108};
   //Double_t data10[nBins]	= {0.864, 0.754, 0.692, 0.568, 0.534, 0.506, 0.428, 0.368, 0.302, 0.242, 0.21, 0.212, 0.204, 0.188, 0.18, 0.116, 0.102};
   //Double_t data11[nBins]	= {0.85, 0.768, 0.678, 0.56, 0.448, 0.45, 0.406, 0.316, 0.288, 0.252, 0.214, 0.192, 0.178, 0.13, 0.144, 0.082, 0.102};
   //Double_t data12[nBins]	= {0.854, 0.778, 0.704, 0.548, 0.464, 0.332, 0.312, 0.296, 0.212, 0.232, 0.204, 0.172, 0.13, 0.132, 0.126, 0.114, 0.07};
   //Double_t data13[nBins]	= {0.826, 0.712, 0.606, 0.546, 0.432, 0.37, 0.312, 0.296, 0.248, 0.22, 0.17, 0.156, 0.12, 0.086, 0.098, 0.086, 0.062};
   //Double_t data14[nBins]	= {0.862, 0.752, 0.566, 0.516, 0.472, 0.372, 0.328, 0.274, 0.2, 0.176, 0.152, 0.148, 0.128, 0.098, 0.066, 0.07, 0.056};
   //Double_t data15[nBins]	= {0.854, 0.732, 0.608, 0.486, 0.388, 0.32, 0.266, 0.238, 0.216, 0.142, 0.136, 0.114, 0.116, 0.07, 0.074, 0.086, 0.044};
   //Double_t data16[nBins]	= {0.81, 0.718, 0.552, 0.49, 0.382, 0.348, 0.314, 0.25, 0.184, 0.128, 0.132, 0.096, 0.084, 0.068, 0.052, 0.048, 0.036};
   //Double_t data17[nBins]	= {0.846, 0.714, 0.576, 0.498, 0.356, 0.308, 0.238, 0.226, 0.144, 0.18, 0.102, 0.08, 0.076, 0.062, 0.066, 0.048, 0.03};
   
   //probability of success |x| as datasets -> klappt nicht
   //Double_t data1[nBins] 	= {0.888, 0.896, 0.884, 0.882, 0.888, 0.862, 0.898, 0.856, 0.864, 0.864, 0.85, 0.854, 0.826, 0.862, 0.854, 0.81, 0.846};
   //Double_t data2[nBins]	= {0.798, 0.818, 0.85, 0.848, 0.83, 0.768, 0.802, 0.76, 0.78, 0.754, 0.768, 0.778, 0.712, 0.752, 0.732, 0.718, 0.714}; 
   //Double_t data3[nBins]	= {0.734, 0.78, 0.782, 0.736, 0.77, 0.702, 0.754, 0.71, 0.678, 0.692, 0.678, 0.704, 0.606, 0.566, 0.608, 0.552, 0.576};
   //Double_t data4[nBins]	= {0.654, 0.734, 0.726, 0.736, 0.7, 0.716, 0.676, 0.634, 0.632, 0.568, 0.56, 0.548, 0.546, 0.516, 0.486, 0.49, 0.498};
   //Double_t data5[nBins]	= {0.648, 0.712, 0.742, 0.696, 0.648, 0.622, 0.586, 0.566, 0.552, 0.534, 0.448, 0.464, 0.432, 0.472, 0.388, 0.382, 0.356};
   //Double_t data6[nBins]	= {0.558, 0.638, 0.63, 0.616, 0.598, 0.598, 0.528, 0.514, 0.488, 0.506, 0.45, 0.332, 0.37, 0.372, 0.32, 0.348, 0.308};
   //Double_t data7[nBins]	= {0.484, 0.622, 0.586, 0.616, 0.526, 0.524, 0.484, 0.504, 0.406, 0.428, 0.406, 0.312, 0.312, 0.328, 0.266, 0.314, 0.238};
   //Double_t data8[nBins]	= {0.46, 0.536, 0.58, 0.558, 0.498, 0.512, 0.458, 0.406, 0.392, 0.368, 0.316, 0.296, 0.296, 0.274, 0.238, 0.25, 0.226};
   //Double_t data9[nBins]	= {0.458, 0.536, 0.516, 0.546, 0.478, 0.466, 0.424, 0.376, 0.35, 0.302, 0.288, 0.212, 0.248, 0.2, 0.216, 0.184, 0.144};
   //Double_t data10[nBins]	= {0.414, 0.538, 0.502, 0.494, 0.466, 0.446, 0.382, 0.34, 0.308, 0.242, 0.252, 0.232, 0.22, 0.176, 0.142, 0.128, 0.18};
   //Double_t data11[nBins]	= {0.388, 0.48, 0.464, 0.434, 0.446, 0.376, 0.328, 0.322, 0.264, 0.21, 0.214, 0.204, 0.17, 0.152, 0.136, 0.132, 0.102};
   //Double_t data12[nBins]	= {0.324, 0.43, 0.456, 0.422, 0.404, 0.352, 0.304, 0.272, 0.228, 0.212, 0.192, 0.172, 0.156, 0.148, 0.114, 0.096, 0.08};
   //Double_t data13[nBins]	= {0.296, 0.45, 0.39, 0.41, 0.344, 0.332, 0.296, 0.2, 0.24, 0.204, 0.178, 0.13, 0.12, 0.128, 0.116, 0.084, 0.076};
   //Double_t data14[nBins]	= {0.28, 0.348, 0.368, 0.396, 0.336, 0.324, 0.25, 0.23, 0.196, 0.188, 0.13, 0.132, 0.086, 0.098, 0.07, 0.068, 0.062};
   //Double_t data15[nBins]	= {0.24, 0.33, 0.366, 0.332, 0.308, 0.31, 0.244, 0.2, 0.18, 0.18, 0.144, 0.126, 0.098, 0.066, 0.074, 0.052, 0.066};
   //Double_t data16[nBins]	= {0.222, 0.37, 0.376, 0.32, 0.244, 0.246, 0.234, 0.184, 0.164, 0.116, 0.082, 0.114, 0.086, 0.07, 0.086, 0.048, 0.048};
   //Double_t data17[nBins]	= {0.246, 0.298, 0.358, 0.334, 0.256, 0.246, 0.216, 0.168, 0.108, 0.102, 0.102, 0.07, 0.062, 0.056, 0.044, 0.036, 0.03};
   
   //Colors
   Int_t i = 1000;
   TColor *color = new TColor(++i, (Double_t)238/(Double_t)255, (Double_t)238/(Double_t)255, (Double_t)238/(Double_t)255);
   TColor *color = new TColor(++i, (Double_t)221/(Double_t)255, (Double_t)221/(Double_t)255, (Double_t)221/(Double_t)255);
   TColor *color = new TColor(++i, (Double_t)204/(Double_t)255, (Double_t)204/(Double_t)255, (Double_t)204/(Double_t)255);
   TColor *color = new TColor(++i, (Double_t)178/(Double_t)255, (Double_t)178/(Double_t)255, (Double_t)178/(Double_t)255);
   TColor *color = new TColor(++i, (Double_t)153/(Double_t)255, (Double_t)153/(Double_t)255, (Double_t)153/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)128/(Double_t)255, (Double_t)128/(Double_t)255, (Double_t)128/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)102/(Double_t)255, (Double_t)102/(Double_t)255, (Double_t)102/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)51/(Double_t)255, (Double_t)51/(Double_t)255, (Double_t)51/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)28/(Double_t)255, (Double_t)28/(Double_t)255, (Double_t)28/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)17/(Double_t)255, (Double_t)17/(Double_t)255, (Double_t)17/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)102/(Double_t)255, (Double_t)51/(Double_t)255, (Double_t)51/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)153/(Double_t)255, (Double_t)0/(Double_t)255, (Double_t)0/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)204/(Double_t)255, (Double_t)0/(Double_t)255, (Double_t)0/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)255/(Double_t)255, (Double_t)102/(Double_t)255, (Double_t)102/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)255/(Double_t)255, (Double_t)153/(Double_t)255, (Double_t)153/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)204/(Double_t)255, (Double_t)153/(Double_t)255, (Double_t)153/(Double_t)255); 
   TColor *color = new TColor(++i, (Double_t)255/(Double_t)255, (Double_t)204/(Double_t)255, (Double_t)204/(Double_t)255); 
   
   //Blue Colors for estimations
   TColor *color = new TColor(1040, (Double_t)114/(Double_t)255, (Double_t)159/(Double_t)255, (Double_t)207/(Double_t)255); 
   TColor *color = new TColor(1060, (Double_t)52/(Double_t)255, (Double_t)101/(Double_t)255, (Double_t)164/(Double_t)255); 
   TColor *color = new TColor(1081, (Double_t)0/(Double_t)255, (Double_t)69/(Double_t)255, (Double_t)134/(Double_t)255); 
   
   //Histograms
   TH1F *hist1 = new TH1F("hist1","",nBins,2,18);
   hist1->SetMarkerStyle(20);
   hist1->SetMarkerSize(0.4);
   hist1->SetMarkerColor(1001);
   hist1->SetStats(0);
   
   TH1F *hist2 = new TH1F("hist2","",nBins,2,18);
   hist2->SetMarkerStyle(20);
   hist2->SetMarkerSize(0.4);
   hist2->SetMarkerColor(1002);
   hist2->SetStats(0);
   
   TH1F *hist3 = new TH1F("hist3","",nBins,2,18);
   hist3->SetMarkerStyle(20);
   hist3->SetMarkerSize(0.4);
   hist3->SetMarkerColor(1003);
   hist3->SetStats(0);

   TH1F *hist4 = new TH1F("hist4","",nBins,2,18);
   hist4->SetMarkerStyle(20);
   hist4->SetMarkerSize(0.4);
   hist4->SetMarkerColor(1004);
   hist4->SetStats(0);
   
   TH1F *hist5 = new TH1F("hist5","",nBins,2,18);
   hist5->SetMarkerStyle(20);
   hist5->SetMarkerSize(0.4);
   hist5->SetMarkerColor(1005);
   hist5->SetStats(0);

   TH1F *hist6 = new TH1F("hist6","",nBins,2,18);
   hist6->SetMarkerStyle(20);
   hist6->SetMarkerSize(0.4);
   hist6->SetMarkerColor(1006);
   hist6->SetStats(0);

   TH1F *hist7 = new TH1F("hist7","",nBins,2,18);
   hist7->SetMarkerStyle(20);
   hist7->SetMarkerSize(0.4);
   hist7->SetMarkerColor(1007);
   hist7->SetStats(0);

   TH1F *hist8 = new TH1F("hist8","",nBins,2,18);
   hist8->SetMarkerStyle(20);
   hist8->SetMarkerSize(0.4);
   hist8->SetMarkerColor(1008);
   hist8->SetStats(0);

   TH1F *hist9 = new TH1F("hist9","",nBins,2,18);
   hist9->SetMarkerStyle(20);
   hist9->SetMarkerSize(0.4);
   hist9->SetMarkerColor(1009);
   hist9->SetStats(0);

   TH1F *hist10 = new TH1F("hist10","",nBins,2,18);
   hist10->SetMarkerStyle(20);
   hist10->SetMarkerSize(0.4);
   hist10->SetMarkerColor(1010);
   hist10->SetStats(0);

   TH1F *hist11 = new TH1F("hist11","",nBins,2,18);
   hist11->SetMarkerStyle(20);
   hist11->SetMarkerSize(0.4);
   hist11->SetMarkerColor(1011);
   hist11->SetStats(0);

   TH1F *hist12 = new TH1F("hist12","",nBins,2,18);
   hist12->SetMarkerStyle(20);
   hist12->SetMarkerSize(0.4);
   hist12->SetMarkerColor(1012);
   hist12->SetStats(0);

   TH1F *hist13 = new TH1F("hist13","",nBins,2,18);
   hist13->SetMarkerStyle(20);
   hist13->SetMarkerSize(0.4);
   hist13->SetMarkerColor(1013);
   hist13->SetStats(0);

   TH1F *hist14 = new TH1F("hist14","",nBins,2,18);
   hist14->SetMarkerStyle(20);
   hist14->SetMarkerSize(0.4);
   hist14->SetMarkerColor(1014);
   hist14->SetStats(0);

   TH1F *hist15 = new TH1F("hist15","",nBins,2,18);
   hist15->SetMarkerStyle(20);
   hist15->SetMarkerSize(0.4);
   hist15->SetMarkerColor(1015);
   hist15->SetStats(0);

   TH1F *hist16 = new TH1F("hist16","",nBins,2,18);
   hist16->SetMarkerStyle(20);
   hist16->SetMarkerSize(0.4);
   hist16->SetMarkerColor(1016);
   hist16->SetStats(0);

   TH1F *hist17 = new TH1F("hist17","",nBins,2,18);
   hist17->SetMarkerStyle(20);
   hist17->SetMarkerSize(0.4);
   hist17->SetMarkerColor(1017);
   hist17->SetStats(0);
   
   //Fill Histograms  
   for(int i=0; i < nBins;  i++)
   {
		if (data1[i] != 0)
			hist1->SetBinContent(i+1,data1[i]);
		if (data2[i] != 0)
			hist2->SetBinContent(i+1,data2[i]);
		if (data3[i] != 0)
			hist3->SetBinContent(i+1,data3[i]);
		if (data4[i] != 0)
			hist4->SetBinContent(i+1,data4[i]);
		if (data5[i] != 0)
			hist5->SetBinContent(i+1,data5[i]);
		if (data6[i] != 0)
			hist6->SetBinContent(i+1,data6[i]);
		if (data7[i] != 0)
			hist7->SetBinContent(i+1,data7[i]);
		if (data8[i] != 0)
			hist8->SetBinContent(i+1,data8[i]);
		if (data9[i] != 0)
			hist9->SetBinContent(i+1,data9[i]);
		if (data10[i] != 0)
			hist10->SetBinContent(i+1,data10[i]);
		if (data11[i] != 0)
			hist11->SetBinContent(i+1,data11[i]);
		if (data12[i] != 0)
			hist12->SetBinContent(i+1,data12[i]);
		if (data13[i] != 0)
			hist13->SetBinContent(i+1,data13[i]);
		if (data14[i] != 0)
			hist14->SetBinContent(i+1,data14[i]);
		if (data15[i] != 0)
			hist15->SetBinContent(i+1,data15[i]);
		if (data16[i] != 0)
			hist16->SetBinContent(i+1,data16[i]);
		if (data17[i] != 0)
			hist17->SetBinContent(i+1,data17[i]);
   }
   
   //Fitting Function
   TF1 *fit = new TF1("fit",fitFunction,2,30,3);
   fit->SetNpx(500);
   fit->SetLineWidth(4);
   
   // Fitting
   TFitResultPtr fitres;
   Double_t par[18][3];
   Double_t err[18][3];
   
   
   fit->SetLineColor(1017);
   hist17->Fit("fit","V+","p ");
   fit->GetParameters(par[17]);
   err[17][0] = fit->GetParError(0);
   err[17][1] = fit->GetParError(1);
   err[17][2] = fit->GetParError(2);
     
   fit->SetLineColor(1001);
   hist1->Fit("fit","V+","p SAME");
   fit->GetParameters(par[1]);
   err[1][0] = fit->GetParError(0);
   err[1][1] = fit->GetParError(1);
   err[1][2] = fit->GetParError(2);
   
   fit->SetLineColor(1002);
   hist2->Fit("fit","V+","p SAME");
   fit->GetParameters(par[2]);
   err[2][0] = fit->GetParError(0);
   err[2][1] = fit->GetParError(1);
   err[2][2] = fit->GetParError(2);
   
   fit->SetLineColor(1003);
   hist3->Fit("fit","V+","p SAME");
   fit->GetParameters(par[3]);
   err[3][0] = fit->GetParError(0);
   err[3][1] = fit->GetParError(1);
   err[3][2] = fit->GetParError(2);
   
   fit->SetLineColor(1004);
   hist4->Fit("fit","V+","p SAME");
   fit->GetParameters(par[4]);
   err[4][0] = fit->GetParError(0);
   err[4][1] = fit->GetParError(1);
   err[4][2] = fit->GetParError(2);
   
   fit->SetLineColor(1005);
   hist5->Fit("fit","V+","p SAME");
   fit->GetParameters(par[5]);
   err[5][0] = fit->GetParError(0);
   err[5][1] = fit->GetParError(1);
   err[5][2] = fit->GetParError(2);
   
   fit->SetLineColor(1006);
   hist6->Fit("fit","V+","p SAME");
   fit->GetParameters(par[6]);
   err[6][0] = fit->GetParError(0);
   err[6][1] = fit->GetParError(1);
   err[6][2] = fit->GetParError(2);
   
   fit->SetLineColor(1007);
   hist7->Fit("fit","V+","p SAME");
   fit->GetParameters(par[7]);
   err[7][0] = fit->GetParError(0);
   err[7][1] = fit->GetParError(1);
   err[7][2] = fit->GetParError(2);
   
   fit->SetLineColor(1008);
   hist8->Fit("fit","V+","p SAME");
   fit->GetParameters(par[8]);
   err[8][0] = fit->GetParError(0);
   err[8][1] = fit->GetParError(1);
   err[8][2] = fit->GetParError(2);
   
   fit->SetLineColor(1009);
   hist9->Fit("fit","V+","p SAME");
   fit->GetParameters(par[9]);
   err[9][0] = fit->GetParError(0);
   err[9][1] = fit->GetParError(1);
   err[9][2] = fit->GetParError(2);
   
   fit->SetLineColor(1010);
   hist10->Fit("fit","V+","p SAME");
   fit->GetParameters(par[10]);
   err[10][0] = fit->GetParError(0);
   err[10][1] = fit->GetParError(1);
   err[10][2] = fit->GetParError(2);
   
   fit->SetLineColor(1011);
   hist11->Fit("fit","V+","p SAME");
   fit->GetParameters(par[11]);
   err[11][0] = fit->GetParError(0);
   err[11][1] = fit->GetParError(1);
   err[11][2] = fit->GetParError(2);
   
   fit->SetLineColor(1012);
   hist12->Fit("fit","V+","p SAME");
   fit->GetParameters(par[12]);
   err[12][0] = fit->GetParError(0);
   err[12][1] = fit->GetParError(1);
   err[12][2] = fit->GetParError(2);
   
   fit->SetLineColor(1013);
   hist13->Fit("fit","V+","p SAME");
   fit->GetParameters(par[13]);
   err[13][0] = fit->GetParError(0);
   err[13][1] = fit->GetParError(1);
   err[13][2] = fit->GetParError(2);
   
   fit->SetLineColor(1014);
   hist14->Fit("fit","V+","p SAME");
   fit->GetParameters(par[14]);
   err[14][0] = fit->GetParError(0);
   err[14][1] = fit->GetParError(1);
   err[14][2] = fit->GetParError(2);
   
   fit->SetLineColor(1015);
   hist15->Fit("fit","V+","p SAME");
   fit->GetParameters(par[15]);
   err[15][0] = fit->GetParError(0);
   err[15][1] = fit->GetParError(1);
   err[15][2] = fit->GetParError(2);
   
   fit->SetLineColor(1016);
   hist16->Fit("fit","V+","p SAME");
   fit->GetParameters(par[16]);
   err[16][0] = fit->GetParError(0);
   err[16][1] = fit->GetParError(1);
   err[16][2] = fit->GetParError(2);

   
   //Draw estimated functions
   TF1 *N40 = new TF1("N40","-1+0.37579*x*x",0,30);
   TF1 *N60 = new TF1("N60","-1+0.54269*x*x",0,30);
   TF1 *N81 = new TF1("N81","-1+0.717935*x*x",0,30);
   
   N40->SetLineColor(1040);
   N60->SetLineColor(1060);
   N81->SetLineColor(1081);
   
   N81->Draw("SAME");
   N60->Draw("SAME");
   N40->Draw("SAME");
   
   for (Int_t i=1; i<18; i++) {
		printf("For $N=%d$ the fitting gives the function: $%f+%f\\cdot x^2$\\\\\n", i+3, par[i][0],par[i][2]);
		//printf("For $|x|=%d$ the fitting gives the function: $%f+%f\\cdot x^2$\\\\\n", i+1, par[i][0],par[i][2]);
		printf("with errors $%f, %f$\\\\\n", err[i][0], err[i][2]);
   }
   printf("\n as csv: \n");
   for (Int_t i=1; i<18; i++) {
		printf("N=%d;%f;%f;%f\n", i+3, par[i][0], par[i][1], par[i][2]);
		//printf("|x|=%d;%f;%f;%f\n", i+1, par[i][0], par[i][1], par[i][2]);
   }
   printf("\n errors \n");	
   for (Int_t i=1; i<18; i++) {	
		printf("N=%d;%f;%f;%f\n", i+3, err[i][0], err[i][1], err[i][2]);
		//printf("|x|=%d;%f;%f;%f\n", i+1, err[i][0], err[i][1], err[i][2]);
   }
   //legend addet in gimp*/
}
