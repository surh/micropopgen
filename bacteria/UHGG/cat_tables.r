# (C) Copyright 2020 Sur Herrera Paredes
# 
# This file is part of micropopgen.
# 
# micropopgen is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# micropopgen is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with micropopgen.  If not, see <http://www.gnu.org/licenses/>.

library(tidyverse)

opts <- commandArgs(trailingOnly = TRUE)
output <- opts[1]
files <- opts[-1]

res <- files %>%
  map_dfr(~read_tsv(.x,
                    col_types = cols(ref_id = col_character(),
                                     ref_pos = col_number(),
                                     coding = col_number(),
                                     .default = col_character())))
write_tsv(res, output)


