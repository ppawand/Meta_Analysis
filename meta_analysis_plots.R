library(ggpubr)


setwd("~/Desktop/Meta-analysis/Plots/")

overalll <-
  ggarrange(overall_biomass_plot,
            overall_activity_plot,
            align = "hv")

ggsave("overall_effects.pdf", overalll, width = 8, height = 4,
       units = "in")

# 
# ecosystem_plots <- 
#   ggarrange(mbc_ecosystem_plot,
#             tbiomass_ecosystem_plot,
#             bacteria_ecosystem_plot,
#             Gpos_ecosystem_plot,
#             Gneg_ecosystem_plot,
#             fungi_ecosystem_plot,
#             tresp_ecosystem_plot,
#             mresp_ecosystem_plot,
#             align = "hv",
#             labels = "auto"
#             )
  

  
mbc_plots <-
  ggarrange(mbc_ecosystem_plot,
           mbc_magnitude_plot,
           mbc_duration_plot,
           mbc_technique_plot,
           mbc_ph_plot,
           align = "hv",
           labels = "auto")

ggsave("mbc_plots.pdf", mbc_plots, width = 12, height = 6,
       units = "in")


bacteria_plots <- 
  ggarrange(bacteria_ecosystem_plot,
            bacteria_magnitude_plot,
            bacteria_duration_plot,
            bacteria_technique_plot,
            bacteria_ph_plot,
            align = "hv",
            labels = "auto")

ggsave("bacteria_plots.pdf", bacteria_plots, width = 12, height = 6,
       units = "in")

fungi_plots <- 
  ggarrange(fungi_ecosystem_plot,
            fungi_magnitude_plot,
            fungi_duration_plot,
            fungi_technique_plot,
            fungi_ph_plot,
            align = "hv",
            labels = "auto")

ggsave("fungi_plots.pdf", fungi_plots, width = 12, height = 6,
       units = "in")

GPositive_plots <-
  ggarrange(Gpos_ecosystem_plot,
            Gpos_magnitude_plot,
            Gpos_duration_plot,
            Gpos_technique_plot,
            Gpos_ph_plot,
            align = "hv",
            labels = "auto")

ggsave("Gpos_bacteria_plots.pdf", GPositive_plots, width = 12, height = 6,
       units = "in")


GNegative_plots <-
  ggarrange(Gneg_ecosystem_plot,
            Gneg_magnitude_plot,
            Gneg_duration_plot,
            Gneg_technique_plot,
            Gneg_ph_plot,
            align = "hv",
            labels = "auto")

ggsave("Gneg_bacteria_plots.pdf", GNegative_plots, width = 12, height = 6,
       units = "in")

TotalBiomass_plots <-
  ggarrange(tbiomass_ecosystem_plot,
            tbiomass_magnitude_plot,
            tbiomass_duration_plot,
            tbiomass_technique_plot,
            tbiomass_ph_plot,
            align = "hv",
            labels = "auto")

ggsave("total_biomass.pdf", TotalBiomass_plots, width = 12, height = 6,
       units = "in")


TotalRespiration_plots <-
  ggarrange(tresp_ecosystem_plot,
            tresp_magnitude_plot,
            tresp_duration_plot,
            tresp_technique_plot,
            tresp_ph_plot,
            align = "hv",
            labels = "auto")

ggsave("Total_respiration.pdf", TotalRespiration_plots, width = 12, height = 6,
       units = "in")

MicrobialRespiration_plots <-
  ggarrange(mresp_ecosystem_plot,
            mresp_magnitude_plot,
            mresp_duration_plot,
            mresp_technique_plot,
            mresp_ph_plot,
            align = "hv",
            labels = "auto")

ggsave("Microbial_respiration.pdf", MicrobialRespiration_plots, width = 12, height = 6,
       units = "in")
 

