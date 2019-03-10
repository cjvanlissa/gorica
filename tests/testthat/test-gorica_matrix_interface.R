hyp <- list(hyp_mat = list(gorica:::parse_hypothesis(c("progAcademic", "progGeneral", "progVocational", "math",
                                                       "progGeneral___X___math", "progVocational___X___math"),
                                                     new_gor$hypotheses[[1]])$hyp_mat[[1]],
                           gorica:::parse_hypothesis(c("progAcademic", "progGeneral", "progVocational", "math", "progGeneral___X___math", "progVocational___X___math"), new_gor$hypotheses[[2]])$hyp_mat[[1]]),
            n_ec = c(3, 0))

gorica(x = model, hypothesis = hyp)
