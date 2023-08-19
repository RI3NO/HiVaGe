install.packages("car")
library(car)
library(dplyr)


# Hiv davies

hiv_var_methods_davies <- data.frame(
  values = c(8.585795737, 66.31516331, 31.0250766, 7.293689327, 26.62542392, 13.72422117, 7.02200577, 9.698939716),
  groups = rep(c("var"))
)

hiv_dist_methods_davies <- data.frame(
  values = c(9.19117766,
             11.60773824,
             6.215706416,
             3.798961976,
             19.5793053,
             11.10345057,
             6.979367037,
             8.424093979,
             10.41005424,
             12.2596289),
  groups = rep(c("dist"))
)

data <- rbind(hiv_var_methods_davies, hiv_dist_methods_davies)

shapiro.test(data["values"][data['groups'] == "var"])

shapiro.test(data["values"][data['groups'] == "dist"])

# Assuming your data is in a data frame called 'data' with a numeric variable 'values' and a categorical variable 'groups'
levene_test_result <- leveneTest(values ~ groups, data = data)
levene_test_result

result <- kruskal.test(values ~ groups, data = data)


var_dist_methods_comparison <- data.frame(dataset="HIV",
                                             metric="davies",
                                             p=result$p.value) 


# Buettner davies

buettner_var_methods_davies <- data.frame(
  values = c(6.33704,
             2.95341,
             2.35916,
             2.55218,
             3.75341,
             2.99143,
             2.86798,
             5.58036,
             6.56752,
             2.76378,
             3.15830),
  groups = rep(c("var"))
)

buettner_dist_methods_davies <- data.frame(
  values = c(6.48255,
             4.53193,
             5.13891,
             6.09004,
             5.82666,
             5.90512,
             2.80018,
             2.16194),
  groups = rep(c("dist"))
)

data <- rbind(buettner_var_methods_davies, buettner_dist_methods_davies)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="buettner",
                                          metric="davies",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)

# Mair davies

mair_var_methods_davies <- data.frame(
  values = c(12.66375203,
             3.284536349,
             5.741518411,
             5.288635653,
             13.00689883,
             5.794134018,
             6.432542406,
             2.790525458,
             4.246771518,
             11.56683722,
             3.321422061),
  groups = rep(c("var"))
)

mair_dist_methods_davies <- data.frame(
  values = c(6.117950249,
             5.076477961,
             5.016136322,
             6.090964724,
             2.354266092,
             5.610900509,
             2.71409589,
             2.403482792,
             5.090777178),
  groups = rep(c("dist"))
)

data <- rbind(mair_var_methods_davies, mair_dist_methods_davies)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="Mair",
                                          metric="davies",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Campbell davies

campbell_var_methods_davies <- data.frame(
  values = c(5.038883041,
             4.174431426,
             4.097662972,
             6.296240832,
             4.881066272,
             5.503309343,
             7.529243778,
             27.23541928),
  groups = rep(c("var"))
)

campbell_dist_methods_davies <- data.frame(
  values = c(3.274148175,
             3.290851589,
             4.020716692,
             3.481717099,
             6.476570945,
             4.340591081,
             4.193663604,
             12.39517418,
             4.067096132),
  groups = rep(c("dist"))
)

data <- rbind(campbell_var_methods_davies, campbell_dist_methods_davies)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="campbell",
                                          metric="davies",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Richard davies

richard_var_methods_davies <- data.frame(
  values = c(3.86069,
             4.77598,
             9.13277,
             4.56138,
             2.55662,
             5.06245,
             5.49795,
             8.57617,
             4.89438,
             4.78169),
  groups = rep(c("var"))
)

richard_dist_methods_davies <- data.frame(
  values = c(6.06986,
             7.27711,
             4.58206,
             4.50913,
             6.84546,
             11.23540,
             7.41128,
             4.16789,
             8.53366,
             10.73675),
  groups = rep(c("dist"))
)

data <- rbind(richard_var_methods_davies, richard_dist_methods_davies)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- t.test(data["values"][data['groups'] == "var"], data["values"][data['groups'] == "dist"])

var_dist_methods_comparison <- data.frame(dataset="richard",
                                          metric="davies",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# All tests conclude that there's not enough strong evidence to assume that variance-based and distribution-based methods differ based on davies


# Calinski
# Hiv 

hiv_var_methods_calinski <- data.frame(
  values = c(8.059862366,
             28.52132615,
             133.6452493,
             91.72878185,
             35.59447354,
             39.60345406,
             163.1434991,
             8.80418086),
  groups = rep(c("var"))
)

hiv_dist_methods_calinski <- data.frame(
  values = c(111.4901625,
             8.485774281,
             48.66412537,
             174.5246756,
             70.93620763,
             7.241442614,
             14.73268027,
             11.41379463,
             16.6034492,
             2.0665228),
  groups = rep(c("dist"))
)

data <- rbind(hiv_var_methods_calinski, hiv_dist_methods_calinski)

shapiro.test(data["values"][data['groups'] == "var"])

shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="HIV",
                                          metric="calinski",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Buettner 

buettner_var_methods_calinski <- data.frame(
  values = c(1.66169,
             23.98358,
             20.38822,
             30.35376,
             13.97297,
             12.04680,
             14.85115,
             2.18394,
             1.03404,
             14.39278,
             11.83871),
  groups = rep(c("var"))
)

buettner_dist_methods_calinski <- data.frame(
  values = c(1.98581,
             2.58552,
             1.93638,
             1.53686,
             1.81206,
             1.58911,
             10.83967,
             17.58947),
  groups = rep(c("dist"))
)

data <- rbind(buettner_var_methods_calinski, buettner_dist_methods_calinski)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="buettner",
                                          metric="calinski",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Mair

mair_var_methods_calinski <- data.frame(
  values = c(22.04735271,
             248.2854217,
             159.4022365,
             168.9121961,
             174.5317843,
             510.1028816,
             570.5494373,
             557.2334778,
             556.8358324,
             32.35380755,
             514.2151361),
  groups = rep(c("var"))
)

mair_dist_methods_calinski <- data.frame(
  values = c(167.46533,
             154.5299786,
             183.4457836,
             151.6819612,
             524.2526783,
             144.8624599,
             565.2004862,
             608.2860232,
             197.8187488),
  groups = rep(c("dist"))
)

data <- rbind(mair_var_methods_calinski, mair_dist_methods_calinski)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="mair",
                                          metric="calinski",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Campbell davies

campbell_var_methods_calinski <- data.frame(
  values = c(593.3815807,
             521.3062344,
             688.9875352,
             524.9778063,
             593.730313,
             696.4782642,
             108.1980569,
             4.385613204),
  groups = rep(c("var"))
)

campbell_dist_methods_calinski <- data.frame(
  values = c(546.7901634,
             807.0706222,
             608.0084106,
             654.1709993,
             556.9836907,
             568.7807633,
             702.8791332,
             143.868893,
             572.5860928),
  groups = rep(c("dist"))
)

data <- rbind(campbell_var_methods_calinski, campbell_dist_methods_calinski)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="campbell",
                                          metric="calinski",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Richard

richard_var_methods_calinski <- data.frame(
  values = c(19.35555,
             15.16258,
             19.43200,
             18.06737,
             7.72145,
             83.71918,
             67.41216,
             11.12494,
             45.00641,
             13.75087),
  groups = rep(c("var"))
)

richard_dist_methods_calinski <- data.frame(
  values = c(15.76576,
             3.18172,
             16.21027,
             15.94631,
             15.24076,
             1.41547,
             11.93063,
             16.03701,
             16.29040,
             1.91250),
  groups = rep(c("dist"))
)

data <- rbind(richard_var_methods_calinski, richard_dist_methods_calinski)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="richard",
                                          metric="calinski",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)



# Average silhouette width
# Hiv 

hiv_var_methods_silhouette <- data.frame(
  values = c(-0.270384145,
             -0.156194235,
             -0.15443785,
             -0.350142112,
             -0.316431571,
             -0.155111916,
             -0.080495508,
             -0.300639636),
  groups = rep(c("var"))
)

hiv_dist_methods_silhouette <- data.frame(
  values = c(-0.140066189,
             -0.007207533,
             -0.299442424,
             0.021934872,
             -0.438387314,
             -0.253294879,
             -0.063025467,
             -0.082104205,
             -0.2669033,
             -0.161916014),
  groups = rep(c("dist"))
)

data <- rbind(hiv_var_methods_silhouette, hiv_dist_methods_silhouette)

shapiro.test(data["values"][data['groups'] == "var"])

shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="HIV",
                                          metric="Average silhouette width",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Buettner 

buettner_var_methods_silhouette <- data.frame(
  values = c(-0.03565,
             0.11048,
             0.07374,
             0.04099,
             0.06842,
             0.06233,
             0.07272,
             0.00290,
             -0.03403,
             0.05474,
             0.08249
  ),
  groups = rep(c("var"))
)

buettner_dist_methods_silhouette <- data.frame(
  values = c(0.00888,
             -0.00356,
             -0.04403,
             0.00104,
             -0.03191,
             -0.00076,
             0.10328,
             0.11500),
  groups = rep(c("dist"))
)

data <- rbind(buettner_var_methods_silhouette, buettner_dist_methods_silhouette)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="buettner",
                                          metric="Average silhouette width",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Mair

mair_var_methods_silhouette <- data.frame(
  values = c(-0.183377706,
             0.047530223,
             -0.041085181,
             -0.023201435,
             -0.10729709,
             -0.296672804,
             -0.077823938,
             -0.034262487,
             0.050951909,
             0.001048907,
             0.025794635),
  groups = rep(c("var"))
)

mair_dist_methods_silhouette <- data.frame(
  values = c(-0.014744744,
             0.001368365,
             -0.020091329,
             -0.024376207,
             0.056152309,
             -0.005616816,
             0.033954603,
             0.03229739,
             0.007554965),
  groups = rep(c("dist"))
)

data <- rbind(mair_var_methods_silhouette, mair_dist_methods_silhouette)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="mair",
                                          metric="Average silhouette width",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Campbell davies

campbell_var_methods_silhouette <- data.frame(
  values = c(-0.104065158,
             -0.102217304,
             -0.083582827,
             -0.308859879,
             -0.109096786,
             -0.126107701,
             -0.329637985,
             -0.527934024),
  groups = rep(c("var"))
)

campbell_dist_methods_silhouette <- data.frame(
  values = c(-0.104226381,
             -0.116680456,
             -0.10329479,
             -0.08464071,
             -0.153468083,
             -0.111211652,
             -0.079592699,
             -0.175677115,
             -0.113189948),
  groups = rep(c("dist"))
)

data <- rbind(campbell_var_methods_silhouette, campbell_dist_methods_silhouette)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="campbell",
                                          metric="Average silhouette width",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Richard

richard_var_methods_silhouette <- data.frame(
  values = c(-0.19869,
             -0.29381,
             -0.14749,
             -0.07482,
             -0.12699,
             -0.11819,
             -0.12739,
             -0.21308,
             -0.13029,
             -0.08859,
             -0.22737),
  groups = rep(c("var"))
)

richard_dist_methods_silhouette <- data.frame(
  values = c(-0.14589,
             -0.19134,
             -0.15012,
             -0.08917,
             -0.06727,
             -0.25651,
             -0.17030,
             -0.27271,
             -0.20914),
  groups = rep(c("dist"))
)

data <- rbind(richard_var_methods_silhouette, richard_dist_methods_silhouette)

shapiro.test(data["values"][data['groups'] == "var"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(dataset="richard",
                                          metric="ilhouette width",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)



setwd("C:/Users/hulum/OneDrive/Documents/GitHub/HiVaGe/Benchmarking/Scripts")

write.csv(var_dist_methods_comparison, file = "var_dist_methods_comparison.csv", row.names = FALSE)




# Average silhouette width

var_methods_silhouette <- data.frame(
  values = c(-0.270384145,
             -0.156194235,
             -0.15443785,
             -0.350142112,
             -0.316431571,
             -0.155111916,
             -0.080495508,
             -0.300639636,
             -0.03565,
             0.11048,
             0.07374,
             0.04099,
             0.06842,
             0.06233,
             0.07272,
             0.00290,
             -0.03403,
             0.05474,
             0.08249,
             -0.183377706,
             0.047530223,
             -0.041085181,
             -0.023201435,
             -0.10729709,
             -0.296672804,
             -0.077823938,
             -0.034262487,
             0.050951909,
             0.001048907,
             0.025794635,
             -0.104065158,
             -0.102217304,
             -0.083582827,
             -0.308859879,
             -0.109096786,
             -0.126107701,
             -0.329637985,
             -0.527934024,
             -0.19869,
             -0.29381,
             -0.14749,
             -0.07482,
             -0.12699,
             -0.11819,
             -0.12739,
             -0.21308,
             -0.13029,
             -0.08859,
             -0.22737),
  groups = rep(c("var"))
)

dist_methods_silhouette <- data.frame(
  values = c(-0.140066189,
             -0.007207533,
             -0.299442424,
             0.021934872,
             -0.438387314,
             -0.253294879,
             -0.063025467,
             -0.082104205,
             -0.2669033,
             -0.161916014,
             0.00888,
             -0.00356,
             -0.04403,
             0.00104,
             -0.03191,
             -0.00076,
             0.10328,
             0.11500,
             -0.014744744,
             0.001368365,
             -0.020091329,
             -0.024376207,
             0.056152309,
             -0.005616816,
             0.033954603,
             0.03229739,
             0.007554965,
             -0.104226381,
             -0.116680456,
             -0.10329479,
             -0.08464071,
             -0.153468083,
             -0.111211652,
             -0.079592699,
             -0.175677115,
             -0.113189948,
             -0.14589,
             -0.19134,
             -0.15012,
             -0.08917,
             -0.06727,
             -0.25651,
             -0.17030,
             -0.27271,
             -0.20914),
  groups = rep(c("dist"))
)

data <- rbind(dist_methods_silhouette, var_methods_silhouette)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- t.test(data["values"][data['groups'] == "var"], data["values"][data['groups'] == "dist"])

var_dist_methods_comparison <- data.frame(metric="Average silhouette width",
                                          p=result$p.value) 

# Adjusted Rand index

var_methods_rand <- data.frame(
  values = c(0.008813358,
             0.011330482,
             0.024307194,
             0.013364952,
             0.024513971,
             0.026613658,
             0.005818407,
             0.000676842,
             0.00181,
             0.00750,
             0.00346,
             0.00193,
             0.02615,
             0.01242,
             0.00038,
             0.01187,
             0.00226,
             0.00625,
             0.00617,
             0.17234007,
             0.478477437,
             0.495511265,
             0.362052995,
             0.16737237,
             0.212924677,
             0.342541297,
             0.399424416,
             0.414314509,
             0.351310228,
             0.395280306,
             0.166692473,
             0.286649161,
             0.238234943,
             0.253536367,
             0.252470027,
             0.238700286,
             0.322676736,
             0.003815267,
             0.21244,
             0.19951,
             0.21679,
             0.00140,
             0.20766,
             0.23028,
             0.21524,
             0.17876,
             0.00312,
             0.21564,
             0.18975),
  groups = rep(c("var"))
)

dist_methods_rand <- data.frame(
  values = c(0.018935337,
             0.030652438,
             0.010456776,
             0.046909382,
             0.013492472,
             0.00611659,
             0.000358504,
             0.027378691,
             0.020061575,
             0.006057324,
             0.00382,
             0.03020,
             0.00508,
             0.01017,
             0.00972,
             0.03685,
             0.00325,
             0.01109,
             0.431663115,
             0.463675008,
             0.376701966,
             0.400010438,
             0.473978792,
             0.46011006,
             0.398412672,
             0.506418605,
             0.437083693,
             0.309433869,
             0.306526991,
             0.22775544,
             0.258286056,
             0.133363494,
             0.244959888,
             0.331421874,
             0.173959951,
             0.306715214,
             0.13212,
             0.23462,
             0.12238,
             0.20322,
             0.00761,
             0.08309,
             0.20197,
             0.14122,
             0.03279
  ),
  groups = rep(c("dist"))
)

data <- rbind(dist_methods_rand, var_methods_rand)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(metric="Adjusted Rand index",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Calinski-Harabasz index

var_methods_calinski <- data.frame(
  values = c(8.059862366,
             28.52132615,
             133.6452493,
             91.72878185,
             35.59447354,
             39.60345406,
             163.1434991,
             8.80418086,
             1.66169,
             23.98358,
             20.38822,
             30.35376,
             13.97297,
             12.04680,
             14.85115,
             2.18394,
             1.03404,
             14.39278,
             11.83871,
             22.04735271,
             248.2854217,
             159.4022365,
             168.9121961,
             174.5317843,
             510.1028816,
             570.5494373,
             557.2334778,
             556.8358324,
             32.35380755,
             514.2151361,
             593.3815807,
             521.3062344,
             688.9875352,
             524.9778063,
             593.730313,
             696.4782642,
             108.1980569,
             4.385613204,
             19.35555,
             15.16258,
             19.43200,
             18.06737,
             7.72145,
             83.71918,
             67.41216,
             11.12494,
             45.00641,
             13.75087),
  groups = rep(c("var"))
)

dist_methods_calinski <- data.frame(
  values = c(111.4901625,
             8.485774281,
             48.66412537,
             174.5246756,
             70.93620763,
             7.241442614,
             14.73268027,
             11.41379463,
             16.6034492,
             2.0665228,
             1.98581,
             2.58552,
             1.93638,
             1.53686,
             1.81206,
             1.58911,
             10.83967,
             17.58947,
             167.46533,
             154.5299786,
             183.4457836,
             151.6819612,
             524.2526783,
             144.8624599,
             565.2004862,
             608.2860232,
             197.8187488,
             546.7901634,
             807.0706222,
             608.0084106,
             654.1709993,
             556.9836907,
             568.7807633,
             702.8791332,
             143.868893,
             572.5860928,
             15.76576,
             3.18172,
             16.21027,
             15.94631,
             15.24076,
             1.41547,
             11.93063,
             16.03701,
             16.29040,
             1.91250),
  groups = rep(c("dist"))
)

data <- rbind(dist_methods_rand, var_methods_rand)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(metric="Calinski-Harabasz index",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)


# Davies-Bouldin index

var_methods_davies <- data.frame(
  values = c(8.585795737,
             66.31516331,
             31.0250766,
             7.293689327,
             26.62542392,
             13.72422117,
             17.02200577,
             9.698939716,
             6.33704,
             2.95341,
             2.35916,
             2.55218,
             3.75341,
             2.99143,
             2.86798,
             5.58036,
             6.56752,
             2.76378,
             3.15830,
             6.33704,
             2.95341,
             2.35916,
             2.55218,
             3.75341,
             2.99143,
             2.86798,
             5.58036,
             6.56752,
             2.76378,
             3.15830,
             5.038883041,
             4.174431426,
             4.097662972,
             6.296240832,
             4.881066272,
             5.503309343,
             7.529243778,
             27.23541928,
             3.86069,
             4.77598,
             9.13277,
             4.56138,
             2.55662,
             5.06245,
             5.49795,
             8.57617,
             4.89438,
             4.78169),
  groups = rep(c("var"))
)

dist_methods_davies <- data.frame(
  values = c(9.19117766,
             11.60773824,
             6.215706416,
             3.798961976,
             19.5793053,
             11.10345057,
             6.979367037,
             8.424093979,
             10.41005424,
             12.2596289,
             6.48255,
             4.53193,
             5.13891,
             6.09004,
             5.82666,
             5.90512,
             2.80018,
             2.16194,
             6.117950249,
             5.076477961,
             5.016136322,
             6.090964724,
             2.354266092,
             5.610900509,
             2.71409589,
             2.403482792,
             5.090777178,
             3.274148175,
             3.290851589,
             4.020716692,
             3.481717099,
             6.476570945,
             4.340591081,
             4.193663604,
             12.39517418,
             4.067096132,
             6.06986,
             7.27711,
             4.58206,
             4.50913,
             6.84546,
             11.23540,
             7.41128,
             4.16789,
             8.53366,
             10.73675),
  groups = rep(c("dist"))
)

data <- rbind(dist_methods_davies, var_methods_davies)

shapiro.test(data["values"][data['groups'] == "var"])
shapiro.test(data["values"][data['groups'] == "dist"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

var_dist_methods_comparison <- data.frame(metric="Davies-Bouldin index",
                                          p=result$p.value) %>%
  rbind(var_dist_methods_comparison,.)

write.csv(var_dist_methods_comparison, file = "var_dist_methods_comparison.csv", row.names = FALSE)




# SIEVE and respective methods

# Average silhouette width

reg_methods_silhouette <- data.frame(
  values = c(-0.2669033,
             -0.063025467,
             -0.253294879,
             -0.15443785,
             0.11048,
             -0.03191,
             0.047530223,
             -0.005616816,
             -0.041085181,
             -0.104065158,
             -0.111211652,
             -0.102217304,
             -0.25651,
             -0.27271,
             -0.19869,
             -0.29381),
  groups = rep(c("regular"))
)

SIEVE_methods_silhouette <- data.frame(
  values = c(0.021934872,
             -0.438387314,
             -0.350142112,
             -0.316431571,
             -0.04403,
             0.04099,
             -0.024376207,
             -0.023201435,
             -0.10729709,
             -0.08464071,
             -0.083582827,
             -0.308859879,
             -0.13029,
             -0.14749,
             -0.07482,
             -0.08917),
  groups = rep(c("SIEVE"))
)

data <- rbind(SIEVE_methods_silhouette, reg_methods_silhouette)

shapiro.test(data["values"][data['groups'] == "regular"])
shapiro.test(data["values"][data['groups'] == "SIEVE"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

reg_SIEVE_methods_comparison <- data.frame(metric="Average silhouette width",
                                          p=result$p.value) 

# Adjusted Rand index

reg_methods_rand <- data.frame(
  values = c(0.020061575,
             0.000358504,
             0.00611659,
             0.024307194,
             0.00750,
             0.00972,
             0.478477437,
             0.46011006,
             0.495511265,
             0.166692473,
             0.244959888,
             0.286649161,
             0.08309,
             0.14122,
             0.21244,
             0.19951),
  groups = rep(c("regular"))
)

SIEVE_methods_rand <- data.frame(
  values = c(0.046909382,
             0.013492472,
             0.013364952,
             0.024513971,
             0.00508,
             0.00193,
             0.400010438,
             0.362052995,
             0.16737237,
             0.258286056,
             0.238234943,
             0.253536367,
             0.00312,
             0.21679,
             0.00140,
             0.20322),
  groups = rep(c("SIEVE"))
)

data <- rbind(SIEVE_methods_rand, reg_methods_rand)

shapiro.test(data["values"][data['groups'] == "regular"])
shapiro.test(data["values"][data['groups'] == "SIEVE"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)
result <- kruskal.test(values ~ groups, data = data)
reg_SIEVE_methods_comparison <- data.frame(metric="Adjusted Rand index",
                                          p=result$p.value) %>%
  rbind(reg_SIEVE_methods_comparison,.)


# Calinski-Harabasz index

reg_methods_calinski <- data.frame(
  values = c(16.6034492,
             14.73268027,
             7.241442614,
             133.6452493,
             23.98358,
             1.81206,
             248.2854217,
             144.8624599,
             159.4022365,
             593.3815807,
             568.7807633,
             521.3062344,
             11.93063,
             16.29040,
             15.16258,
             19.43200),
  groups = rep(c("regular"))
)

SIEVE_methods_calinski <- data.frame(
  values = c(48.66412537,
             8.485774281,
             35.59447354,
             39.60345406,
             2.58552,
             30.35376,
             154.5299786,
             168.9121961,
             174.5317843,
             807.0706222,
             688.9875352,
             524.9778063,
             3.18172,
             18.06737,
             7.72145,
             16.21027),
  groups = rep(c("SIEVE"))
)

data <- rbind(SIEVE_methods_rand, reg_methods_rand)

shapiro.test(data["values"][data['groups'] == "regular"])
shapiro.test(data["values"][data['groups'] == "SIEVE"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

reg_SIEVE_methods_comparison <- data.frame(metric="Calinski-Harabasz index",
                                          p=result$p.value) %>%
  rbind(reg_SIEVE_methods_comparison,.)


# Davies-Bouldin index

reg_methods_davies <- data.frame(
  values = c(10.41005424,
             6.979367037,
             11.10345057,
             31.0250766,
             2.95341,
             5.82666,
             3.284536349,
             5.610900509,
             5.741518411,
             5.038883041,
             4.340591081,
             4.174431426,
             7.41128,
             8.53366,
             4.77598,
             9.13277),
  groups = rep(c("regular"))
)

SIEVE_methods_davies <- data.frame(
  values = c(6.215706416,
             11.60773824,
             26.62542392,
             13.72422117,
             4.53193,
             2.55218,
             5.076477961,
             5.288635653,
             13.00689883,
             3.290851589,
             4.097662972,
             6.296240832,
             7.27711,
             4.56138,
             2.55662,
             4.58206),
  groups = rep(c("SIEVE"))
)

data <- rbind(SIEVE_methods_davies, reg_methods_davies)

shapiro.test(data["values"][data['groups'] == "regular"])
shapiro.test(data["values"][data['groups'] == "SIEVE"])

leveneTest(values ~ groups, data = data)

result <- kruskal.test(values ~ groups, data = data)

reg_SIEVE_methods_comparison <- data.frame(metric="Davies-Bouldin index",
                                          p=result$p.value) %>%
  rbind(reg_SIEVE_methods_comparison,.)

write.csv(reg_SIEVE_methods_comparison, file = "reg_SIEVE_methods_comparison.csv", row.names = FALSE)


# t-test of purity of scVEGs

scVEGs_purity <- 0.469387755

all_purities <- c(0.314176245,
                  0.30573711,
                  0.313705234,
                  0.32124031,
                  0.328372093,
                  0.31874612,
                  0.310454401,
                  0.32744186,
                  0.331782946,
                  0.320524836,
                  0.332093023,
                  0.354166667,
                  0.330062152,
                  0.320500165,
                  0.299566295,
                  0.325898389,
                  0.321251549)

shapiro.test(all_purities)

result <- t.test(all_purities, mu = scVEGs_purity, alternative = "two.sided")

scVEGs_comparison <- data.frame(metric="Purity", p=result$p.value) 

# test of Calinski-Harabasz index of scVEGs

scVEGs_calinski <- 91.72878185

all_calinski <- c(16.6034492,
                  2.0665228,
                  8.059862366,
                  14.73268027,
                  11.41379463,
                  28.52132615,
                  7.241442614,
                  133.6452493,
                  174.5246756,
                  70.93620763,
                  48.66412537,
                  8.485774281,
                  35.59447354,
                  39.60345406,
                  163.1434991,
                  8.80418086,
                  111.4901625)

shapiro.test(all_calinski)

library(e1071)
library(moments)
library(tse)
install.packages('tseries')

jarque.bera.test(all_calinski)

result <- wilcox.test(all_calinski, mu = scVEGs_calinski, alternative = "two.sided")

scVEGs_comparison <- data.frame(metric="Calinski-Harabasz index", p=result$p.value) %>%
  rbind(scVEGs_comparison,.)

# test of Davies-Bouldin index index of scVEGs

scVEGs_davies <- 7.293689327

all_davies <- c(10.41005424,
                12.2596289,
                8.585795737,
                6.979367037,
                8.424093979,
                66.31516331,
                11.10345057,
                31.0250766,
                3.798961976,
                19.5793053,
                6.215706416,
                11.60773824,
                26.62542392,
                13.72422117,
                17.02200577,
                9.698939716,
                9.19117766)

shapiro.test(all_davies)

jarque.bera.test(all_davies)

# There is no test for it

# test of Average silhouette width index of scVEGs

scVEGs_silhouette <- -0.300639636

all_silhouettes <- c(-0.2669033,
                -0.161916014,
                -0.270384145,
                -0.063025467,
                -0.082104205,
                -0.156194235,
                -0.253294879,
                -0.15443785,
                -0.007207533,
                -0.299442424,
                0.021934872,
                -0.438387314,
                -0.350142112,
                -0.316431571,
                -0.155111916,
                -0.080495508,
                -0.140066189)

shapiro.test(all_silhouettes)

result <- t.test(all_silhouettes, mu = scVEGs_silhouette, alternative = "two.sided")

scVEGs_comparison <- data.frame(metric="Average silhouette width index", p=result$p.value) %>%
  rbind(scVEGs_comparison,.)

# test of Average rand index index of scVEGs

scVEGs_rand <- 0.000676842

all_rand <- c(0.020061575,
                     0.006057324,
                     0.008813358,
                     0.000358504,
                     0.027378691,
                     0.011330482,
                     0.00611659,
                     0.024307194,
                     0.030652438,
                     0.010456776,
                     0.046909382,
                     0.013492472,
                     0.013364952,
                     0.024513971,
                     0.026613658,
                     0.005818407,
                     0.018935337)

shapiro.test(all_rand)

result <- t.test(all_rand, mu = scVEGs_rand, alternative = "two.sided")

scVEGs_comparison <- data.frame(metric="Average rand index", p=result$p.value) %>%
  rbind(scVEGs_comparison,.)

write.csv(scVEGs_comparison, file = "scVEGs_comparison.csv", row.names = FALSE)


# M3Drop test on Richard

# test of purity

M3Drop_purity <- 0.56124

all_purities <- c(0.39583,
                  0.57828,
                  0.58081,
                  0.56566,
                  0.38413,
                  0.57513,
                  0.57513,
                  0.57008,
                  0.57828,
                  0.53914,
                  0.57386,
                  0.80702,
                  0.56629,
                  0.57895,
                  0.57323,
                  0.57891,
                  0.57955,
                  0.57702,
                  0.58018)

shapiro.test(all_purities)

jarque.bera.test(all_purities)

#  No test as data's not symmetrical and not normally distributed

# test of Calinski-Harabasz index of M3Drop

M3Drop_calinski <- 16.29040

all_calinski <- c(1.91250,
                  11.93063,
                  16.03701,
                  15.16258,
                  1.41547,
                  19.43200,
                  15.94631,
                  16.21027,
                  15.76576,
                  11.12494,
                  15.24076,
                  3.18172,
                  18.06737,
                  7.72145,
                  19.35555,
                  45.00641,
                  13.75087,
                  83.71918,
                  67.41216)

shapiro.test(all_calinski)

jarque.bera.test(all_calinski)

#  No test as data's not symmetrical and not normally distributed


# test of Davies-Bouldin index of M3Drop

M3Drop_davies <- 8.53366

all_davies <- c(10.73675,
                7.41128,
                4.16789,
                4.77598,
                11.23540,
                9.13277,
                4.50913,
                4.58206,
                6.06986,
                8.57617,
                6.84546,
                7.27711,
                4.56138,
                2.55662,
                3.86069,
                4.89438,
                4.78169,
                5.06245,
                5.49795)

shapiro.test(all_davies)

result <- t.test(all_davies, mu = M3Drop_davies, alternative = "less")

M3Drop_comparison <- data.frame(metric="Davies-Bouldin index", p=result$p.value)


# test of Average silhouette width index of M3Drop

M3Drop_silhouette <- -0.27271

all_silhouettes <- c(-0.20914,
                     -0.25651,
                     -0.17030,
                     -0.19869,
                     -0.06727,
                     -0.29381,
                     -0.19134,
                     -0.08917,
                     -0.14589,
                     -0.21308,
                     -0.15012,
                     -0.13029
                     -0.14749,
                     -0.07482,
                     -0.08859,
                     -0.22737,
                     -0.12699,
                     -0.11819,
                     -0.12739)

shapiro.test(all_silhouettes)

result <- t.test(all_silhouettes, mu = M3Drop_silhouette, alternative = "less")

M3Drop_comparison <- data.frame(metric="Average silhouette width index", p=result$p.value) %>%
  rbind(M3Drop_comparison,.)

# test of Average rand index index of M3Drop

M3Drop_rand <- 0.14122

all_rand <- c(0.03279,
              0.08309,
              0.20197,
              0.21244,
              0.00761,
              0.19951,
              0.23462,
              0.20322,
              0.13212,
              0.17876,
              0.12238,
              0.00312,
              0.21679,
              0.00140,
              0.21564,
              0.18975,
              0.20766,
              0.23028,
              0.21524)

shapiro.test(all_rand)

jarque.bera.test(all_rand)

result <- wilcox.test(all_rand, mu = M3Drop_rand, alternative = "less")

M3Drop_comparison <- data.frame(metric="Average rand index", p=result$p.value) %>%
  rbind(M3Drop_comparison,.)

write.csv(M3Drop_comparison, file = "M3Drop_comparison.csv", row.names = FALSE)



# Methods ranking

norm_df <- data.frame(matrix(nrow = 22, ncol = 0))

install.packages("readxl")
library(readxl)

data <- read_excel("statistical tests.xlsx", sheet = "Sheet5")

data <- as.data.frame(data)

rownames(data) <- data[,1]

data <- data[,2:33]

data[,1:26]

class(normalized_scores)
# All metrix excepd davies 
for (i in 1:26) {
  data_column <- na.omit(data[,i])
  mean_score <- mean(data_column)
  sd_score <- sd(data_column)
  normalized_scores <- (data_column - mean_score) / sd_score
  normalized_scores <- as.matrix(normalized_scores)
  norm_df <- rbind(as.data.frame(normalized_scores), matrix(NA, nrow =(22-length(normalized_scores)), ncol = 1)) %>%
    cbind(norm_df,.)
}


data[,27:32]

for (i in 27:32) {
  data_column <- na.omit(data[,i])
  mean_score <- mean(data_column)
  sd_score <- sd(data_column)
  normalized_scores <- (data_column - mean_score) / sd_score
  normalized_scores <- as.matrix(-normalized_scores)
  norm_df <- rbind(as.data.frame(normalized_scores), matrix(NA, nrow =(22-length(normalized_scores)), ncol = 1)) %>%
    cbind(norm_df,.)
}

write.csv(norm_df, "z_norm_values.csv", row.names = FALSE)

data <- read_excel("z_norm_values.xlsx", sheet = "z_norm_values")

data <- as.data.frame(data)

rownames(data) <- data[,1]

data <- data[,2:33]

Average_metrics <- as.data.frame(matrix(nrow = 22, ncol = 1))
rownames(Average_metrics) <- rownames(data)
for (i in 1:32) {
  Average_metrics[i,] <- mean(na.omit(as.numeric(data[i,])))
  }


-Average_metrics[order(Average_metrics$V1), ]

write.csv(Average_metrics, "Average_metrics.csv", row.names = TRUE)

for (i in 1:22) {
  filename <- paste("group", i, sep = "")
  print(filename)
}


class(data[2,])

group1 <- na.omit(data[1,])
group2 <- na.omit(data[2,])
group3 <- na.omit(data[3,])
group4 <- na.omit(data[4,])
group5 <- na.omit(data[5,])
group6 <- na.omit(data[6,])
group7 <- na.omit(data[7,])
group8 <- na.omit(data[8,])
group9 <- na.omit(data[9,])
group10 <- na.omit(data[10,])
group11 <- na.omit(data[11,])
group12 <- na.omit(data[12,])
group13 <- na.omit(data[13,])
group14 <- na.omit(data[14,])
group15 <- na.omit(data[15,])
group16 <- na.omit(data[16,])
group17 <- na.omit(data[17,])
group18 <- na.omit(data[18,])
group19 <- na.omit(data[19,])
group20 <- na.omit(data[20,])
group21 <- na.omit(data[21,])
group22 <- na.omit(data[22,])

rownames(data)
data_column <- na.omit(data[,i])

group_data <- list(group22, group9, group10, group5, group12, group18, group17, group8, group6, group13, group11, group21, group1, group2, group19, group4, group15, group7, group14, group16, group20, group3)
                   
par(las = 2)  
boxplot(transposed_df, col = "#c5dcf2", names = rep("", length(transposed_df[1,])))

transposed_df <- t(data)


column_means <- c()

for (i in 1:22) {
  column_means <- c(column_means, mean(na.omit(transposed_df[,i]) )  
}

# Get the indices of columns sorted by means in descending order
sorted_indices <- order(column_means, decreasing = TRUE)

rearranged_df <- transposed_df[, sorted_indices]

boxplot(rearranged_df, col = "#c5dcf2", names = rep("", length(transposed_df[1,])))

