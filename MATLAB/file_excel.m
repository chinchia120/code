score = xlsread("score_data.xlsx");
score = xlsread("score_data.xlsx", "B2: D4");

M = mean(score')';
xlswrite("score_data.xlsx", M, 1, "E2: E4") %file_name, variable, sheet, range
xlswrite("score_data.xlsx", {'Mean'}, 1, "E1");

[score header] = xlsread("score_data.xlsx") %[numeric string]