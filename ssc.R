#### setup ####
require(tidyverse)


#### load data as db ####
# db <-

#### Venn diagram ####
ggVennDiagram::ggVennDiagram(x = list(
  Lungs     = which(db$Polmone == 1),
  Heart     = which(db$Cuore == 1),
  Kidneys   = which(db$Rene == 1),
  Esophagus = which(db$Esofago == 1)
)) + guides(fill = FALSE)

#### modeling ####
model_esof <- glm(Esofago ~ g_età + g_disdur + g_ttt + g_cute_rodnan +
                    `a_anti-SSA`+`a_anti-RNP`+`a_anti-Rna polimerasi III`+
                    `a_anti-Scl 70`+`a_anti-Centromero`+`a_Pm-Scl 100/75`,
                  data = db, family = "binomial")

model_cuor <- glm(Cuore ~ g_età + g_disdur + g_ttt + g_cute_rodnan + 
                    `a_anti-SSA`+`a_anti-RNP`+`a_anti-Rna polimerasi III`+
                    `a_anti-Scl 70`+`a_anti-Centromero`+`a_Pm-Scl 100/75`,
                  data = db, family = "binomial")

model_polm <- glm(Polmone ~ g_età + g_disdur + g_ttt + g_cute_rodnan + 
                    `a_anti-SSA`+`a_anti-RNP`+`a_anti-Rna polimerasi III`+
                    `a_anti-Scl 70`+`a_anti-Centromero`+`a_Pm-Scl 100/75`,
                  data = db, family = "binomial")

model_rene <- glm(Rene ~ g_età + g_disdur + g_ttt + g_cute_rodnan + 
                    `a_anti-SSA`+`a_anti-RNP`+`a_anti-Rna polimerasi III`+
                    `a_anti-Scl 70`+`a_anti-Centromero`+`a_Pm-Scl 100/75`,
                  data = db, family = "binomial")

#### predictions ####
predictions_horiz <- select(db, outcomes)
predictions_horiz$pred_esof <- predict(model_esof, type = "response")
predictions_horiz$pred_cuor <- predict(model_cuor, type = "response")
predictions_horiz$pred_polm <- predict(model_polm, type = "response")
predictions_horiz$pred_rene <- predict(model_rene, type = "response")

predictions <- predictions_horiz %>%
  rowid_to_column("id") %>%
  unite("Esofago", Esofago, pred_esof) %>%
  unite("Cuore",   Cuore,   pred_cuor) %>%
  unite("Polmone", Polmone, pred_polm) %>%
  unite("Rene",    Rene,    pred_rene) %>%
  gather("site", "outcome", outcomes) %>%
  separate(outcome, c("outcome", "prediction"), sep = "_", convert = TRUE)


#### ROC curves ####
roc_esof   <- pROC::roc(predictions_horiz, Esofago, pred_esof)
ggroc_esof <- pROC::ggroc(roc_esof) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 1, lty = 2) +
  ggtitle(paste("Esophagus, AUC =", round(pROC::auc(roc_esof), 3)))

roc_polm   <- pROC::roc(predictions_horiz, Polmone, pred_polm)
ggroc_polm <- pROC::ggroc(roc_polm) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 1, lty = 2) +
  ggtitle(paste("Lungs, AUC =", round(pROC::auc(roc_polm), 3)))

roc_rene   <- pROC::roc(predictions_horiz, Rene, pred_rene)
ggroc_rene <- pROC::ggroc(roc_rene) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 1, lty = 2) +
  ggtitle(paste("Kidney, AUC =", round(pROC::auc(roc_rene), 3)))

roc_cuor   <- pROC::roc(predictions_horiz, Cuore, pred_cuor)
ggroc_cour <- pROC::ggroc(roc_cuor) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 1, lty = 2) +
  ggtitle(paste("Heart, AUC =", round(pROC::auc(roc_cuor), 3)))

four_rocs <- ggpubr::ggarrange(ggroc_esof, ggroc_cour, ggroc_polm, ggroc_rene, ncol = 2, nrow = 2)
# ggsave(four_rocs, "roc.pdf", width = 11)

roc_all  <- pROC::roc(predictions, outcome, prediction)
ggroc_all <- pROC::ggroc(roc_all) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 1, lty = 2) +
  ggtitle(paste("Overall, AUC =", round(pROC::auc(roc_all), 3)))

ggpubr::ggarrange(four_rocs, ggroc_all)


#### Goodness of fit ####
predictions %>%
  group_by(id) %>%
  arrange(prediction) %>%
  mutate(score = row_number() * outcome) %>%
  mutate(score = ifelse(score > 0, score - 1, score)) %>% 
  mutate(max_score = case_when(
    sum(outcome) == 1 ~ 3,
    sum(outcome) == 2 ~ 5,
    TRUE ~ 6)) %>%
  summarise(score = sum(score), max_score = first(max_score)) %>%
  mutate(final = score/max_score) %>%
  summarise(performance = mean(final))

# one row for each predicted complication
high_risk_db <- predictions %>%
  group_by(id) %>%
  arrange(desc(prediction)) %>% 
  mutate(n_tot = sum(outcome)) %>%
  filter(row_number() <= n_tot) %>%
  select(-outcome, risk_group = site) %>%
  left_join(rowid_to_column(db, "id"))

caret::confusionMatrix(
  data = factor(high_risk_db_vert$prediction),
  reference = factor(high_risk_db_vert$outcome))
