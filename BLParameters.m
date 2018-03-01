function [ P, Omega, q ] = BLParameters(tau, mktNoShares, periodReturns, currentPrices)
    NoViews = 6;
    P = zeros(NoViews, size(mktNoShares,2)); 
    Omega = zeros(size(mktNoShares,2), size(mktNoShares,2)); 
    q = zeros(NoViews,1);
    Q = cov(periodReturns);

    % 1st View: Pertaining to Company side. Here we assumge that our 10 largest
    % market cap stocks will underperform our 10 smallest cap stocks by 2%

    % We need to first sort (and keep track of indices) for our stocks
    [order, I] = sort(mktNoShares.*currentPrices);

    for i=1:size(order,1)
        if I(i) > 10
            P(1,i) = 1;
        else
            P(1,i) = -1;
        end
    end

    % no idea why I need to do this but I apparently do
    clear oder clear I

    % Outperformance of our ten largest by 3%
    q(1) = 0.03;

    % We calculate our confidence in this view
    Omega(1,1) = tau * P(1,:) * Q * P(1,:)'; 

    % 2nd View: Pertaining to the sensitivity in the relationship between the 
    % 'Information Technology' and 'Financials Sector'. If the combined average 
    % return of AAPL and IBM exceeds the combined average return of C, WFC and 
    % JPM, our view will be 25% higher than the difference between the average. 
    % If underperformed, the view will be 10% lower.
    P(2,8:10) = -ones(1,3);
    P(2,11:12) = ones(1,2);

    avg_performance = mean(geomean(periodReturns(:,8:10)+1)-1) - mean(geomean(periodReturns(:,11:12)+1)-1);

    if avg_performance >= 0
        amplifcation = 1.25;
    else
        amplifcation = 1.1;
    end

    q(2) = amplifcation * avg_performance;

    % Confidence Level
    Omega(2,2) = tau * P(2,:) * Q * P(2,:)';

    % 3rd and 4th Views: Pertaining to momentum in Energy stocks.  If the average
    % return of XOM and MRO is positive, the view predicts the average return
    % will be 5% higher. If negative, 15% lower.
    P(3,15) = 1;
    q(3) = 1.05*(geomean(periodReturns(:,15)+1)-1);
    Omega(3,3) = tau * P(3,:) * Q * P(3,:)';

    P(4,16) = 1;
    q(4) = 0.85*(geomean(periodReturns(:,16)+1)-1);
    Omega(4,4) = tau * P(4,:) * Q * P(4,:)';

    % 5th and 6th Views: Pertaining to momentum in Telecommunications stocks. 
    % If the average return of T and Vz is positive, the view predicts the 
    % average return will be 10% higher. If negative, 1% lower.
    P(5,19) = 1;
    q(5) = 1.1*(geomean(periodReturns(:,19)+1)-1);
    Omega(5,5) = tau * P(5,:) * Q * P(5,:)';

    P(6,20) = 1;
    q(6) = 0.99*(geomean(periodReturns(:,20)+1)-1);
    Omega(6,6) = tau * P(6,:) * Q * P(6,:)';

end

