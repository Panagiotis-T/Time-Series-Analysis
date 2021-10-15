🏎️ Formula 1 Lap Times Prediction

This project models Formula 1 drivers’ lap times in a particular race using statistical methods, accounting for driver-specific patterns, pit stops, and variable lap-time variance.

📊 Overview

An initial OLS model predicts the parameters and forecasts lap times, including uncertainty.

For some drivers, the GLM model shows large autocorrelations in residuals, systematically underestimating lap times before pit stops and overestimating afterward.

Even after correcting for spikes, lap-time variance remains higher during specific laps.

To address this, WLS (Weighted Least Squares) estimates are applied, factoring in correlations between lap times and increased variance during events causing spikes.

The model parameters are re-estimated by assuming larger standard deviation during pit stop laps.

The optimal covariance matrix is determined using a relaxation algorithm, estimating correlation and three parameters related to lap-time spikes: first laps, laps with pit stops, and laps following pit stops.

🔹 ARMA & Seasonal Processes

The project also explores ARMA processes and seasonal effects, demonstrating how the choice of coefficients in operator polynomials affects the structure of the process.

Simulated data and empirical autocorrelation functions are used for illustration.

📁 Data

The data, along with data from other races, is publically available at: https: //www.kaggle.com/rohanrao/formula-1-world-championship-1950-2020/version/13?select=lap_times.csv.
