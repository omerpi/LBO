clc; clear; close all;

%% Inputs
L = 0;

% Step 1: Assumptions
years_w_data = 6;
pur_yr = 2023;
exit_yr = 2028;
mul.pur_ev_ebitda = 7.25;
mul.exit_ev_ebitda = 7;
pur_ebitda = 1714;
curr_debt = -1982;
cash_on_balance_sheet = 1902;

duration = exit_yr - pur_yr; 

% Step 2: Sources & Uses
% Sources
sources.loans.revolver.cap = 750;
sources.loans.revolver.ebitda_mul = 0;
sources.loans.revolver.int_rate = 0.065*ones(1, duration);

sources.loans.term_a.ebitda_mul = 1.5;
sources.loans.term_a.int_rate = (L + 0.02)*ones(1, duration);

sources.loans.term_b.ebitda_mul = 2.2;
sources.loans.term_b.int_rate = (L + 0.03)*ones(1, duration);

sources.loans.senior_notes.ebitda_mul = 2;
sources.loans.senior_notes.int_rate = 0.095*ones(1, duration);

sources.loans.junior_notes.ebitda_mul = 0;
sources.loans.junior_notes.int_rate = 0.09*ones(1, duration);

sources.loans.other_debt.ebitda_mul = 0;
sources.loans.other_debt.int_rate = 0.11*ones(1, duration);

sources.equity.managment_equity = 360;

% Uses
transaction_fees_and_expenses_prc = 0.02;
cash_for_ops = 1300;
other_expenses = 0;


% Step 3a: Free Cash Flow (past data)
revenue = [8590.5 8436.6 8750.7 9584 12293.4 12368.2];
adj_ebitda_past = [769.5 684.2 760.0 1212.9 2382.2 1843.1];
capex_past = [-476.7 -198.2 -218.8 -224.2 -353.7 -407.2];
net_wc_change_past = [101.2 111.4 -200.1 707.0 -221.3 -545.7];


% Step 3b: Free Cash Flow less Interest Expense before Debt repayment
DnA = [-237.7 -243.8 -270.4 -290.0	-322.6 -365.5];

revenue_growth_future = 0.04*ones(1, duration);
ebitda_margin_future = 0.15*ones(1, duration);
DnA_prc_of_sales_future = -0.03*ones(1, duration);
capex_prc_of_sales_future = -0.032*ones(1, duration);
delta_nwc_prc_of_sales_future = -0.1*ones(1, duration);
tax_prc.future_est = -0.23*ones(1, duration);


% Step 5: Debt Balances
SOFR.future_est = 0*ones(1, duration);
SOFR.floor = 0.054;

sources.loans.revolver.cash_sweep = 1;
sources.loans.term_a.cash_sweep = 1;
sources.loans.term_b.cash_sweep = 1;
sources.loans.senior_notes.cash_sweep = 1;
sources.loans.junior_notes.cash_sweep = 1;
sources.loans.other_debt.cash_sweep = 1;

sources.loans.revolver.amortization = 0*ones(1, duration);
sources.loans.term_a.amortization = 0.2*ones(1, duration);
sources.loans.term_b.amortization = 0.01*ones(1, duration);
sources.loans.senior_notes.amortization = 0*ones(1, duration);
sources.loans.junior_notes.amortization = 0*ones(1, duration);
sources.loans.other_debt.amortization = 0*ones(1, duration);


%% STEP 1: Transaction Assumptions
pur_ev = pur_ebitda * mul.pur_ev_ebitda;
equity_pur_price = pur_ev + curr_debt + cash_on_balance_sheet;

%% Step 3a: Free Cash Flow
for i = 1:duration
    revenue = [revenue revenue(end)*(1 + revenue_growth_future(i))];
end
adj_ebitda_future = revenue(years_w_data+1:end).*ebitda_margin_future;%%%
adj_ebitda = [adj_ebitda_past adj_ebitda_future];
ebitda_margin = adj_ebitda./revenue;


%% Step 3b: Free Cash Flow less Interest Expense before Debt repayment
revenue_growth = 0;%%%
for i = 2:duration + years_w_data
    revenue_growth = [revenue_growth (revenue(i)-revenue(i-1))/revenue(i-1)];
end

delta_nwc_prc_of_sales = 0;
for i = 2:years_w_data
    delta_nwc_prc_of_sales = [delta_nwc_prc_of_sales net_wc_change_past(i)/(revenue(i)-revenue(i-1))];
end
delta_nwc_prc_of_sales = [delta_nwc_prc_of_sales delta_nwc_prc_of_sales_future];

net_wc_change = net_wc_change_past;
for i = 1:duration
    net_wc_change = [net_wc_change ...
        (revenue(years_w_data+i)-revenue(years_w_data+i-1))*delta_nwc_prc_of_sales(years_w_data+i)];
end

%%%%%%%%% STEP 1
exit_ebidta = adj_ebitda(end); %%%
exit_ev = exit_ebidta * mul.exit_ev_ebitda;

capex_prc_of_sales = [capex_past./revenue(1:years_w_data) capex_prc_of_sales_future];
capex = capex_prc_of_sales.*revenue;

DnA_prc_of_sales = [DnA./revenue(1:years_w_data) DnA_prc_of_sales_future];

DnA = [DnA revenue(years_w_data+1:end).*DnA_prc_of_sales(years_w_data+1:end)];
%% STEP 2: Sources & Uses
% Uses
transaction_fees_and_expenses = transaction_fees_and_expenses_prc * pur_ev;

uses.pur_price = equity_pur_price ;
uses.refinance_existing_debt = -curr_debt;
uses.transaction_fees_and_expenses = transaction_fees_and_expenses;
uses.cash_for_ops = cash_for_ops;
uses.other = other_expenses;
uses.total = 0;
uses.total = sum(struct2array(uses));

% Sources
sources.loans.revolver.amount = sources.loans.revolver.cap*sources.loans.revolver.ebitda_mul;
sources.loans.term_a.amount = pur_ebitda*sources.loans.term_a.ebitda_mul;
sources.loans.term_b.amount = pur_ebitda*sources.loans.term_b.ebitda_mul;
sources.loans.senior_notes.amount = pur_ebitda*sources.loans.senior_notes.ebitda_mul;
sources.loans.junior_notes.amount = pur_ebitda*sources.loans.junior_notes.ebitda_mul;
sources.loans.other_debt.amount = pur_ebitda*sources.loans.other_debt.ebitda_mul;

sources.debt_ebitda_mul = 0;
sources.total_loans = 0;
loans_fn = fieldnames(sources.loans);
for k = 1:numel(loans_fn)
    curr_field = loans_fn{k};
    if isfield(sources.loans.(curr_field), 'ebitda_mul')
        sources.debt_ebitda_mul = sources.debt_ebitda_mul + sources.loans.(curr_field).ebitda_mul;
        sources.total_loans = sources.total_loans + sources.loans.(curr_field).amount;
    end
end
sources.total = uses.total;
sources.excess_cash = cash_on_balance_sheet;
sources.equity.sponser_equity = sources.total - sources.total_loans - sources.equity.managment_equity...
    - sources.excess_cash;



%% Back to FCF:
SOFR.avg_val = max((SOFR.future_est - [0 SOFR.future_est(2:end)])/2, SOFR.floor);


% Balances After Mandatory
loans_fn_partial = loans_fn;
revolver_index = find(strcmp(loans_fn_partial, 'revolver'));
loans_fn_partial(revolver_index) = [];

for k = 1:numel(loans_fn_partial)
    curr_field = loans_fn_partial{k};
    debt_balances.after_mandatory_repay.(curr_field) = sources.loans.(curr_field).amount;
    
    for jj = 1:duration
        curr_amortizaion_val = sources.loans.(curr_field).amount * sources.loans.(curr_field).amortization(jj);
        debt_balances.after_mandatory_repay.(curr_field) = [debt_balances.after_mandatory_repay.(curr_field) ...
            max(debt_balances.after_mandatory_repay.(curr_field)(end) - curr_amortizaion_val, 0)];
        sources.loans.(curr_field).pmt = curr_amortizaion_val;
    end
end

%%
debt_balances.total_debt = zeros(1, duration+1);
for k = 1:numel(loans_fn)
    curr_field = loans_fn{k};
    debt_balances.after_optional_repay.(curr_field) = zeros(1, duration+1);
    debt_balances.after_optional_repay.(curr_field)(1) = sources.loans.(curr_field).amount;
end

revolver = [sources.loans.revolver.amount zeros(1,duration)];

debt_after_mndt = [0*revolver;
    debt_balances.after_mandatory_repay.term_a;...
    debt_balances.after_mandatory_repay.term_b;...
    debt_balances.after_mandatory_repay.senior_notes;...
    debt_balances.after_mandatory_repay.junior_notes;...
    debt_balances.after_mandatory_repay.other_debt];

debt_after_optn = [revolver; debt_after_mndt(2:end, :)];

int_rates = [sources.loans.revolver.int_rate;...
    sources.loans.term_a.int_rate + SOFR.avg_val;...
    sources.loans.term_b.int_rate + SOFR.avg_val;...
    sources.loans.senior_notes.int_rate;...
    sources.loans.junior_notes.int_rate;...
    sources.loans.other_debt.int_rate];

%%
% Calc CFADS
int_expense = zeros(1, duration);
excess_fcf = zeros(1, duration);
pre_tax_pl = zeros(1, duration);
taxes = zeros(1, duration);
fcf_before_debt = zeros(1, duration);

for ii = 1:duration
    if ii == 1
        num_iter = 1;
    else
        num_iter = 10;
    end
    for iter = 1:num_iter
        temp_excess_fcf = 1;
        diff = 0;
        int_expense(ii) = -sum((debt_after_optn(:,ii)/2 + debt_after_optn(:,ii+1)/2).*int_rates(:,ii));
        if ii==1 && 1  %patch due to error in excel
            int_expense(ii) = int_expense(ii) - (debt_after_mndt(2,ii)/2 + debt_after_mndt(2,ii+1)/2)*int_rates(1,ii)...
                - (debt_after_mndt(1,ii)/2 + debt_after_mndt(1,ii+1)/2).*int_rates(1,ii);
        end
    
        pre_tax_pl(ii) = adj_ebitda(years_w_data+ii) + int_expense(ii) + DnA(years_w_data+ii);
        taxes(ii) = pre_tax_pl(ii)*tax_prc.future_est(ii);
        fcf_before_debt(ii) = adj_ebitda(years_w_data+ii) + capex(years_w_data+ii) + net_wc_change(years_w_data+ii) + int_expense(ii) + taxes(ii);
    
        excess_fcf(ii) = fcf_before_debt(ii) + sum(debt_after_mndt(2:end,ii+1)) - sum(debt_after_optn(2:end,ii));
        revolver(ii + 1) = min(sources.loans.revolver.cap, max(0,revolver(ii) - sources.loans.revolver.cash_sweep*excess_fcf(ii)));
        debt_after_mndt(1, ii + 1) = revolver(ii + 1);
        debt_after_optn(1, ii + 1) = revolver(ii + 1);
    
        %Term A
        if sources.loans.term_a.cash_sweep == 1 && excess_fcf(ii) > 0 && revolver(ii) - revolver(ii+1) > -1e-5 && debt_after_optn(2, ii) - sources.loans.term_a.pmt > 1e-5
            debt_after_optn(2, ii + 1) = max(0, debt_after_mndt(2, ii + 1) - excess_fcf(ii) + revolver(ii) - revolver(ii+1));
            if debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1) > 1e-5
                temp_excess_fcf = - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
                     + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1);
            end
            excess_fcf(ii) = fcf_before_debt(ii) + sum(debt_after_mndt(2:end,ii+1)) - sum(debt_after_optn(2:end,ii));
            revolver(ii + 1) = min(sources.loans.revolver.cap, max(0,revolver(ii) - sources.loans.revolver.cash_sweep*excess_fcf(ii)));
            debt_after_mndt(1, ii + 1) = revolver(ii + 1);
            debt_after_optn(1, ii + 1) = revolver(ii + 1);
            int_expense(ii) = -sum((debt_after_optn(:,ii)/2 + debt_after_optn(:,ii+1)/2).*int_rates(:,ii));
            pre_tax_pl(ii) = adj_ebitda(years_w_data+ii) + int_expense(ii) + DnA(years_w_data+ii);
            taxes(ii) = pre_tax_pl(ii)*tax_prc.future_est(ii);
            fcf_before_debt(ii) = adj_ebitda(years_w_data+ii) + capex(years_w_data+ii)+ net_wc_change(years_w_data+ii) + int_expense(ii) + taxes(ii);
    
        else
            debt_after_optn(2, ii + 1) = debt_after_mndt(2, ii + 1);
        end

        %Term B
        if sources.loans.term_b.cash_sweep == 1 && excess_fcf(ii) > 0 && revolver(ii) - revolver(ii+1) > -1e-5 && debt_after_optn(3, ii) > 1e-5 - sources.loans.term_b.pmt && temp_excess_fcf > 1e-5
            debt_after_optn(3, ii + 1) = max(0, debt_after_mndt(3, ii + 1) - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
               + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1));
            if debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1) > 1e-5
                temp_excess_fcf = - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
                    + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1)...
                    + debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1);
            end
            excess_fcf(ii) = fcf_before_debt(ii) + sum(debt_after_mndt(2:end,ii+1)) - sum(debt_after_optn(2:end,ii));
            revolver(ii + 1) = min(sources.loans.revolver.cap, max(0,revolver(ii) - sources.loans.revolver.cash_sweep*excess_fcf(ii)));
            debt_after_mndt(1, ii + 1) = revolver(ii + 1);
            debt_after_optn(1, ii + 1) = revolver(ii + 1);
            int_expense(ii) = -sum((debt_after_optn(:,ii)/2 + debt_after_optn(:,ii+1)/2).*int_rates(:,ii));
            pre_tax_pl(ii) = adj_ebitda(years_w_data+ii) + int_expense(ii) + DnA(years_w_data+ii);
            taxes(ii) = pre_tax_pl(ii)*tax_prc.future_est(ii);
            fcf_before_debt(ii) = adj_ebitda(years_w_data+ii) + capex(years_w_data+ii)+ net_wc_change(years_w_data+ii) + int_expense(ii) + taxes(ii);
    
        else
            debt_after_optn(3, ii + 1) = debt_after_mndt(3, ii + 1);
        end
        %Senior Notes
        if sources.loans.senior_notes.cash_sweep == 1 && excess_fcf(ii) > 0 && revolver(ii) - revolver(ii+1) > -1e-5 && debt_after_optn(4, ii) - sources.loans.senior_notes.pmt > 1e-5 && temp_excess_fcf > 1e-5
            debt_after_optn(4, ii + 1) = max(0, debt_after_mndt(4, ii + 1) - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
               + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1) + debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1));
            if debt_after_mndt(4, ii + 1) - debt_after_optn(4, ii + 1) > 1e-5
                temp_excess_fcf = - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
                    + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1)...
                    + debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1)...
                    + debt_after_mndt(4, ii + 1) - debt_after_optn(4, ii + 1);
            end
            excess_fcf(ii) = fcf_before_debt(ii) + sum(debt_after_mndt(2:end,ii+1)) - sum(debt_after_optn(2:end,ii));
            revolver(ii + 1) = min(sources.loans.revolver.cap, max(0,revolver(ii) - sources.loans.revolver.cash_sweep*excess_fcf(ii)));
            debt_after_mndt(1, ii + 1) = revolver(ii + 1);
            debt_after_optn(1, ii + 1) = revolver(ii + 1);
            int_expense(ii) = -sum((debt_after_optn(:,ii)/2 + debt_after_optn(:,ii+1)/2).*int_rates(:,ii));
            pre_tax_pl(ii) = adj_ebitda(years_w_data+ii) + int_expense(ii) + DnA(years_w_data+ii);
            taxes(ii) = pre_tax_pl(ii)*tax_prc.future_est(ii);
            fcf_before_debt(ii) = adj_ebitda(years_w_data+ii) + capex(years_w_data+ii)+ net_wc_change(years_w_data+ii) + int_expense(ii) + taxes(ii);
        else
            debt_after_optn(4, ii + 1) = debt_after_mndt(4, ii + 1);
        end

        %Junior Notes
        if sources.loans.junior_notes.cash_sweep == 1 && excess_fcf(ii) > 0 && revolver(ii) - revolver(ii+1) > -1e-5 && debt_after_optn(5, ii) - sources.loans.junior_notes.pmt > 1e-5 && temp_excess_fcf > 1e-5
            debt_after_optn(5, ii + 1) = max(0, debt_after_mndt(5, ii + 1) - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
               + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1) + debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1)...
               + debt_after_optn(4, ii + 1) - debt_after_mndt(4, ii + 1));
            if debt_after_mndt(5, ii + 1) - debt_after_optn(5, ii + 1)  > 1e-5
                temp_excess_fcf = - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
                    + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1)...
                    + debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1) ...
                    + debt_after_mndt(4, ii + 1) - debt_after_optn(4, ii + 1) ...
                    + debt_after_mndt(5, ii + 1) - debt_after_optn(5, ii + 1);
                 diff = debt_after_mndt(5, ii + 1) - debt_after_optn(5, ii + 1);
            end
            excess_fcf(ii) = fcf_before_debt(ii) + sum(debt_after_mndt(2:end,ii+1)) - sum(debt_after_optn(2:end,ii));
            revolver(ii + 1) = min(sources.loans.revolver.cap, max(0,revolver(ii) - sources.loans.revolver.cash_sweep*excess_fcf(ii)));
            debt_after_mndt(1, ii + 1) = revolver(ii + 1);
            debt_after_optn(1, ii + 1) = revolver(ii + 1);
            int_expense(ii) = -sum((debt_after_optn(:,ii)/2 + debt_after_optn(:,ii+1)/2).*int_rates(:,ii));
            pre_tax_pl(ii) = adj_ebitda(years_w_data+ii) + int_expense(ii) + DnA(years_w_data+ii);
            taxes(ii) = pre_tax_pl(ii)*tax_prc.future_est(ii);
            fcf_before_debt(ii) = adj_ebitda(years_w_data+ii) + capex(years_w_data+ii)+ net_wc_change(years_w_data+ii) + int_expense(ii) + taxes(ii);
        else
            debt_after_optn(5, ii + 1) = debt_after_mndt(5, ii + 1);
        end

        %Other Debt
        if sources.loans.other_debt.cash_sweep == 1 && excess_fcf(ii) > 0 && revolver(ii) - revolver(ii+1) > -1e-5 && debt_after_optn(6, ii) - sources.loans.other_debt.pmt > 1e-5 && temp_excess_fcf > 1e-5
            debt_after_optn(6, ii + 1) = max(0, debt_after_mndt(6, ii + 1) - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
               + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1) + debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1)...
               + debt_after_optn(4, ii + 1) - debt_after_mndt(4, ii + 1) + debt_after_optn(5, ii + 1) - debt_after_mndt(5, ii + 1));
            if debt_after_mndt(6, ii + 1) - debt_after_optn(6, ii + 1) > 1e-5
                temp_excess_fcf = - excess_fcf(ii) + revolver(ii) - revolver(ii+1)...
                    + debt_after_mndt(2, ii + 1) - debt_after_optn(2, ii + 1)...
                    + debt_after_mndt(3, ii + 1) - debt_after_optn(3, ii + 1)...
                    + debt_after_mndt(4, ii + 1) - debt_after_optn(4, ii + 1)...
                    + debt_after_mndt(5, ii + 1) - debt_after_optn(5, ii + 1)...
                    + debt_after_mndt(6, ii + 1) - debt_after_optn(6, ii + 1);
            end
            excess_fcf(ii) = fcf_before_debt(ii) + sum(debt_after_mndt(2:end,ii+1)) - sum(debt_after_optn(2:end,ii));
            revolver(ii + 1) = min(sources.loans.revolver.cap, max(0,revolver(ii) - sources.loans.revolver.cash_sweep*excess_fcf(ii)));
            debt_after_mndt(1, ii + 1) = revolver(ii + 1);
            debt_after_optn(1, ii + 1) = revolver(ii + 1);
            int_expense(ii) = -sum((debt_after_optn(:,ii)/2 + debt_after_optn(:,ii+1)/2).*int_rates(:,ii));
            pre_tax_pl(ii) = adj_ebitda(years_w_data+ii) + int_expense(ii) + DnA(years_w_data+ii);
            taxes(ii) = pre_tax_pl(ii)*tax_prc.future_est(ii);
            fcf_before_debt(ii) = adj_ebitda(years_w_data+ii) + capex(years_w_data+ii)+ net_wc_change(years_w_data+ii) + int_expense(ii) + taxes(ii);
        else
            debt_after_optn(6, ii + 1) = debt_after_mndt(6, ii + 1);
        end
    end
end


%% Finalizing
total_debt = sum(debt_after_optn);

cash_balance_after_optn_repay = [cash_for_ops zeros(1, duration)];
for ii = 1:duration
    cash_balance_after_optn_repay(ii + 1) = cash_balance_after_optn_repay(ii) + fcf_before_debt(ii)...
    + total_debt(ii + 1) -total_debt(ii);
end

net_debt = total_debt - cash_balance_after_optn_repay;
exit_year_net_debt = net_debt(end);
implied_equity_value = exit_ev - exit_year_net_debt;
management_stake_prc = sources.equity.managment_equity...;
    /(sources.equity.sponser_equity + sources.equity.managment_equity);
management_stake = management_stake_prc*implied_equity_value;
sponser_equity_at_exit = implied_equity_value - management_stake;

sponser_equity_check = sources.equity.sponser_equity;
sponser_equity_to_sources = sources.equity.sponser_equity/sources.total;

% Calc IRR & MOIC
sponser_irr = power((sponser_equity_at_exit/sponser_equity_check),1/duration) - 1;
sponser_moic = sponser_equity_at_exit/sponser_equity_check;

%% Credit Stats
net_debt_to_adj_ebitda = net_debt(2:end)./adj_ebitda(years_w_data + 1:end);
senior_secured_debt_to_adj_ebitda = sum(debt_after_optn(1:3,2:end))./adj_ebitda(years_w_data + 1:end);
adj_ebitda_to_interest = adj_ebitda(years_w_data + 1:end)./(-int_expense);
adj_ebitda_less_capex_to_interest = (adj_ebitda(years_w_data + 1:end) + capex(years_w_data + 1:end))...
    ./(-int_expense);

%% Print
disp("Sponser's IRR: " + sponser_irr*100 + "%");
disp("Sponser's Cash-on-Cash: " + sponser_moic + "X");


