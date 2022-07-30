
if choose_set_figures == 1 % TVP-VARX in the main text
    
    cholesky        = 0;     % Cholesky identification
    istore          = 1;%1;     % Store year-by-year impulse responses
    kappa_Q         = 0.06;%0.015; % Prior: time-variation in B
    
    % Start Training Sample
    start_year      = 1975;
    start_month     = 1;
    % Start TVP-VAR
    start_year_TVP  = 1991;
    start_month_TVP = 1;
    end_year        = 2017;
    end_month       = 9;

    % Choose Period IRF confidence bands
    year_confidence     = 2003;
    month_confidence    = 1;
    % Choose Periods for comparison
    year_compare_stock_1  = 1991;
    month_compare_stock_1 = 1;
    year_compare_stock_2  = 2002;
    month_compare_stock_2 = 1;
    year_compare_stock_3  = 2007;
    month_compare_stock_3 = 12;

    year_compare_house_1  = 1991;
    month_compare_house_1 = 1;
    year_compare_house_2  = 1995;
    month_compare_house_2 = 1;
    year_compare_house_3  = 2007;
    month_compare_house_3 = 12;

    year_compare_sacr_1   = 1995;
    month_compare_sacr_1  = 1;
    year_compare_sacr_2   = 2007;
    month_compare_sacr_2  = 12;

    % Proxy
    instr_use             = 1;
    
elseif choose_set_figures == 2 % TVP-VAR Cholesky identification
    
    cholesky        = 1;     % Cholesky identification
    istore          = 0;     % Store year-by-year impulse responses
    kappa_Q         = 0.01;  % Prior: time-variation in B
    
    % Start Training Sample
    start_year      = 1978;
    start_month     = 11;
    % Start TVP-VAR
    start_year_TVP  = 1985; %1991;
    start_month_TVP = 1;
    end_year        = 2017;
    end_month       = 9;
    
    % Proxy
    instr_use       = 1;
    
elseif choose_set_figures == 3 % TVP-VARX sensitivity "Priors"
    
    cholesky        = 0;     % Cholesky identification
    istore          = 0;     % Store year-by-year impulse responses
    kappa_Q         = 0.01;  % Prior: time-variation in B
    
    % Start Training Sample
    start_year      = 1978;
    start_month     = 11;
    % Start TVP-VAR
    start_year_TVP  = 1985; %1991;
    start_month_TVP = 1;
    end_year        = 2017;
    end_month       = 9;
    
    % Proxy
    instr_use       = 1;
    
elseif choose_set_figures == 4 % TVP-VARX sensitivity "ELB Episode"
    
    cholesky        = 0;     % Cholesky identification
    istore          = 0;     % Store year-by-year impulse responses
    kappa_Q         = 0.015; % Prior: time-variation in B
    
    % Start Training Sample
    start_year      = 1978;
    start_month     = 11;
    % Start TVP-VAR
    start_year_TVP  = 1990;
    start_month_TVP = 1;
    end_year        = 2007;
    end_month       = 12;

    % Proxy
    instr_use       = 1;
    
elseif choose_set_figures == 5 % TVP-VARX sensitivity "Timing of Actions"
    
    cholesky        = 0;     % Cholesky identification
    istore          = 0;     % Store year-by-year impulse responses
    kappa_Q         = 0.01;  % Prior: time-variation in B
    
    % Start Training Sample
    start_year      = 1978;
    start_month     = 11;
    % Start TVP-VAR
    start_year_TVP  = 1991;
    start_month_TVP = 1;
    end_year        = 2017;
    end_month       = 9;
    
    % Proxy
    instr_use       = 2;
    
end

if colored_figures == 1
    irf_color   = 'Mediumblue';
    band_color  = 'blue';
    dash_color  = 'Crimson';
else
    irf_color   = 'Black';
    band_color  = 'k';
    dash_color  = 'DimGrey';
end
