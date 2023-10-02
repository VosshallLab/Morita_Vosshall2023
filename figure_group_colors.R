#assign colors and labels to groups in each experiment type, keep consistent across figures
#run before each script

#orco
orco <- c()
orco$col <- c("black",
              "gray30",
              "gray50",
              "firebrick3")

orco$fill <- c("black",
              "gray30",
              "gray50",
              "firebrick3")

orco$breaks = c("orl_wt",
                "orco_5_het",
                "orco_16_het",
                "orco_516_ko")

orco$labels = c("Wildtype",
                expression(italic("orco" ^ "5/+")),
                expression(italic("orco" ^ "16/+")),
                expression(italic(paste("orco"^"5/16"))))

#Ir25a
ir25a <- c()
ir25a$col <- c("black",
               "gray30",
               "gray50",
               "dodgerblue3")

ir25a$fill <- c("black",
               "gray30",
               "gray50",
               "dodgerblue3")

ir25a$breaks = c("lvp_wt",
                 "ir25a_BamHI_het",
                 "ir25a_19_het",
                 "ir25a_BamHI19_ko")

ir25a$labels = c("Wildtype",
                 expression(italic("Ir25a" ^ "BamHI/+")),
                 expression(italic("Ir25a" ^ "19/+")),
                 expression(italic(paste("Ir25a"^"BamHI/19"))))

#Ir76b
ir76b <- c()
ir76b$col <- c("black",
               "gray30",
               "gray50",
               "royalblue4")

ir76b$fill <- c("black",
               "gray30",
               "gray50",
               "royalblue4")


ir76b$breaks = c("lvp_wt",
                 "ir76b_61_het",
                 "ir76b_32_het",
                 "ir76b_6132_ko")

ir76b$labels = c("Wildtype",
                 expression(italic("Ir76b" ^ "61/+")),
                 expression(italic("Ir76b" ^ "32/+")),
                 expression(italic(paste("Ir76b"^"61/32"))))

#Ir8a
ir8a <- c()
ir8a$col <- c("black",
               "gray30",
               "gray50",
               "orchid4")

ir8a$fill <- c("black",
              "gray30",
              "gray50",
              "orchid4")

ir8a$breaks = c("orl_wt",
                 "ir8a_dsred_het",
                 "ir8a_attp_het",
                 "ir8a_dsredattp_ko")

ir8a$labels = c("Wildtype",
                 expression(italic("Ir8a" ^ "DsRed/+")),
                 expression(italic("Ir8a" ^ "attp/+")),
                 expression(italic(paste("Ir8a"^"DsRed/attp"))))


#antenna cut
antcut <- c()
antcut$col <- c("gray70",
                "black",
                "black",
                "black")

antcut$fill <- c("gray70",
                "black",
                "gray60",
                "white")

antcut$breaks = c("noheat",
                  "control",
                  "antennatipcut",
                  "antennafullcut")

antcut$labels = c("No heat",
                  "Intact",
                  "Tip removed",
                  "All removed")

#tarsi cut
tarsicut <- c()
tarsicut$col <- c("black",
                  "black",
                  "black",
                  "black")

tarsicut$fill <- c("black",
                  "white",
                  "gray80",
                  "gray80")

tarsicut$breaks = c("CONTROL",
                    "FL",
                    "ML",
                    "HL")

tarsicut$labels = c("Intact",
                    "Foreleg",
                    "Midleg",
                    "Hindleg")

#ir140 orco double mutant
ir140_orco <- c()
ir140_orco$fill <- c("black",
                  "firebrick3",
                  "firebrick3",
                  "firebrick3",
                  "firebrick3")

ir140_orco$col <- c("black",
                   "firebrick3",
                   "gray30",
                   "gray50",
                   "slateblue3")

ir140_orco$breaks = c("orl",
                    "orco1616",
                    "Ir140_14_orco1616",
                    "Ir140_49_orco1616",
                    "Ir140_4914_orco1616")

ir140_orco$labels = c("orl",
                    "orco1616",
                    "Ir140_14_orco1616",
                    "Ir140_49_orco1616",
                    "Ir140_1449_orco1616")

#ir140 single mutant
ir140 <- c()
ir140$col <- c("black",
                    "gray30",
                    "gray50",
                    "slateblue3")

ir140$fill <- c("black",
                     "gray30",
                     "gray50",
                     "slateblue3")

ir140$breaks = c("orl",
                 "Ir140_17_het",
                 "Ir140_144_het",
                 "Ir140_17144_ko")

ir140$labels = c("orl",
                 "Ir140_17_het",
                 "Ir140_144_het",
                 "Ir140_17144_ko")


#combine as one variable
figcolors <- list(orco=orco,
                  ir25a=ir25a,
                  ir76b=ir76b,
                  ir8a=ir8a,
                  antcut=antcut,
                  tarsicut=tarsicut,
                  ir140_orco=ir140_orco,
                  ir140=ir140)

rm(orco,ir25a,ir76b,ir8a,antcut,tarsicut,ir140_orco)
