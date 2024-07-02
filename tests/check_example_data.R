# Test spar ----

## CW data driven random projection ----
check <- function(...) stopifnot(...)
data("example_data")
set.seed(123)
spar_res <- spar(example_data$x,example_data$y,
                 xval=example_data$xtest,yval=example_data$ytest,nummods=c(5,10,15,20,25,30))
coefs <- coef(spar_res)
pred <- predict(spar_res,xnew=example_data$x)

check(all.equal(spar_res$val_res$Meas[1:9],
                c(32013.2551462143, 31987.6169393502, 32019.2383253907,
                  32033.673691955, 31930.5119683019, 32003.1155244146,
                  31616.3540654455, 32022.209440531, 32017.3858154688)))

check(all.equal(coefs$beta[1:9],
                c(0, 0, 0.986224823721168,
                  0, 0, -0.712549431503832,
                  0.336435350851393, 0.264603047060951, 0)))

check(all.equal(pred[1:9], c(-6.10303963576124, -21.7000296025596, -13.3177035327418,
                             13.4827449924309, 16.3988081443907, 52.2097119407757,
                             5.45163765797581, -10.0694269535631, -6.27629028810337)))

plot(spar_res)
plot(spar_res,"Val_Meas","nummod")
plot(spar_res,"Val_numAct","lambda")
plot(spar_res,"coefs",prange=c(1,400))


## CW random projection ----
set.seed(123)
spar_res <- spar(example_data$x,example_data$y,
                 type.rpm = "cw",
                 xval=example_data$xtest,yval=example_data$ytest,nummods=c(5,10,15,20,25,30))
coefs <- coef(spar_res)
pred <- predict(spar_res,xnew=example_data$x)

check(all.equal(spar_res$val_res$Meas[1:9],
                c(32313.2047120677, 32271.4468744033, 32412.5258807132,
                  32531.9042576578, 32423.9058838029, 32386.2599679008,
                  32606.031169539, 32281.6163349935, 31941.4251320425)))

check(all.equal(coefs$beta[1:9],
                c(0, 0, 0.142867443968229, 0,
                  -0.000841312867941294, 0, 0.331078038555172, 0, 0)))

check(all.equal(pred[1:9], c(3.64960972043469, -7.0532526877315, -10.5616523697266,
                             2.27239452978243, 16.9635061270789, 48.9170717937091,
                             3.2406931772349, -22.9414868693306, -13.2665610501716)))


## Achiloptas random projection ----
set.seed(123)
spar_res <- spar(example_data$x,example_data$y,
                 type.rpm = "sparse",
                 xval=example_data$xtest,yval=example_data$ytest,nummods=c(5,10,15,20,25,30))
coefs <- coef(spar_res)
pred <- predict(spar_res,xnew=example_data$x)
check(all.equal(spar_res$val_res$Meas[1:9],
                c(35883.719094479, 35902.7142507063, 35854.7739398024,
                  35893.2699324101, 35633.6246122969, 35612.639051566,
                  35432.0515536489, 35067.3691936891, 34650.0347715329)))

check(all.equal(coefs$beta[1:9],
                c(0.0596462190284467, 0.109222324641079, 0.237057797096033,
                  0.11128478322278, -0.0348209408291419, -0.278058254617255,
                  0.361863864236041, -0.00294441431757295, 0)))

check(all.equal(pred[1:9], c(2.19243426411014, -9.12183174701782, -13.3020470836854,
                             2.01104709113758, 16.047570767589, 49.0194095879441,
                             1.28620453404856, -19.8453887766184, -13.7577326371216)))

## Gaussian random projection ----
set.seed(123)
spar_res <- spar(example_data$x,example_data$y,
                 type.rpm = "gaussian",
                 xval=example_data$xtest,yval=example_data$ytest,nummods=c(5,10,15,20,25,30))
coefs <- coef(spar_res)
pred <- predict(spar_res,xnew=example_data$x)
check(all.equal(spar_res$val_res$Meas[1:9],
                c(31767.9629949645, 31710.3794862918, 31799.471132618,
                  31858.8652255291, 31934.7357733949, 31626.6949852092,
                  31602.0562754849, 32078.3982464915, 32456.7857095489)))

check(all.equal(coefs$beta[1:9],
                c(0.0357910396013952, 0.0398682969230409, 0.397086638535963,
                  0, 0.0461359481400526, -0.345496707960359,
                  0.0906839461630321, 0.372196676857705, 0.0630583767592219)))

check(all.equal(pred[1:9], c(0.955574625383119, -7.52015175441735, -12.6030859260428,
                             1.29554027174273, 14.346591829296, 47.6147673089992,
                             -3.57712863562895, -21.1376344306622, -7.29159259461283)))


# Test spar.cv ----


## CW data driven random projection ----
data("example_data")
set.seed(123)
spar_res <- spar.cv(example_data$x,example_data$y,
                    type.rpm = "cwdatadriven",
                    nummods=c(5,10,15,20,25,30))
check(all.equal(spar_res$betas@x[1:9],
                c(0.052367149047699, 0.00478355893105122, 0.0234420488129777,
                  -0.00172687468161232, -0.0171928987979096, 0.0377331517165617,
                  0.019822890201749, 0.020772310143975, -0.0135123137034562)))

## CW random projection ----
set.seed(123)
spar_res <- spar.cv(example_data$x,example_data$y,nummods=c(5,10,15,20,25,30),
                    type.rpm = "cw")
check(all.equal(spar_res$betas@x[1:9],
                c(-0.0256082521065418, 0.0522571908738069, 0.0196807843841702,
                  0.0459709403122155, -0.0289555937543286, -0.0294509954023284,
                  0.0403494653040291, 0.0155434654476443, 0.0444118781392775)))



## Achiloptas random projection ----
set.seed(123)
spar_res <- spar.cv(example_data$x,example_data$y,nummods=c(5,10,15,20,25,30),
                    type.rpm = "sparse")
check(all.equal(spar_res$betas@x[1:9],
                c(0.0126587829253347, -0.0188845238363487, 0.0399932803339893,
                  0.0109525857427377, -0.000525059406709394, 0.0368063867601782,
                  0.000246764884650729, 0.0277934337912715, -0.0071758939237829)))


## Gaussian random projection ----
set.seed(123)
spar_res <- spar.cv(example_data$x,example_data$y,nummods=c(5,10,15,20,25,30),
                 type.rpm = "gaussian")
check(all.equal(spar_res$betas@x[1:9],
                c(-0.00446497193077751, 0.033361358294459, 0.0561732697026508,
                  -0.0127116665014228, -0.0376739622890854, -0.0310112713849973,
                  0.00262425209145074, 0.0100322291982222, -0.014057846161789)))#

