import numpy as np
from sklearn.linear_model import LinearRegression

import matplotlib.pyplot as plt

# Generate some example data
np.random.seed(0)
X = 2 * np.random.rand(100, 1)
y = 4 + 3 * X + np.random.randn(100, 1)

# Create a linear regression model
model = LinearRegression()
model.fit(X, y)

# Make predictions
X_new = np.array([[0], [2]])
y_predict = model.predict(X_new)

# Plot the data and the regression line
plt.scatter(X, y, color='blue', label='Data points')
plt.plot(X_new, y_predict, color='red', linewidth=2, label='Regression line')
plt.xlabel('X')
plt.ylabel('y')
plt.title('Linear Regression Example')
plt.legend()
plt.show()