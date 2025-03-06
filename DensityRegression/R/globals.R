# Defining global variables for usage in data.tables without getting note
# "no visible binding for global variable..." during R CMD check

utils::globalVariables(c("..covCol", "..densi", "..inds", "..regressors",
                         "..selection", "..single", "..var_vec",
                         "..y", "all_of", "gam_offsets", "gam_weights",
                         "group_id", "ind response", "unweighting_factor",
                         "weighting_factor", ".", "ind", "response", "k",
                         "value", "variable", "..var_ind", "discrete",
                         "..response_name", "Delta", "height", "freq",
                         "weighted_counts", "type"))
