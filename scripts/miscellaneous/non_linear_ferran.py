from lmfit import Minimizer, Parameters
import numpy as np
import matplotlib.pyplot as plt


class NonLinearFitting(object):

    def __init__(self, obj_func, params, x, y):
        self.obj_func = obj_func
        self.params = params
        self.x = x
        self.y = y
        self.out = None

    def least_squares(self, **kwargs):
        """
        Returns:
            res: dict of fitted values
            predicted: array of float
        """
        mini = Minimizer(self.obj_func, self.params, fcn_args=(self.x, self.y))
        self.out = mini.least_squares(self.params)
        predicted = self.y + self.out.residual
        res = {}
        for k in self.out.__dict__['params']:
            res[self.out.params[k].name] = self.out.params[k].value
        return res, predicted

    def minimize(self, **kwargs):
        """
        Returns:
            res: dict of fitted values
            predicted: array of float
        """
        mini = Minimizer(self.obj_func, self.params, fcn_args=(self.x, self.y))
        self.out = mini.minimize(method='cobyla')
        predicted = self.y + self.out.residual
        res = {}
        for k in self.out.__dict__['params']:
            res[self.out.params[k].name] = self.out.params[k].value
        return res, predicted

    def plot_fit(self, **kwargs):
        fitted_curve = self.y + self.out.residual
        plt.plot(self.x, self.y, 'o')
        plt.plot(self.x, fitted_curve, 'r')
        plt.xlabel('x')
        plt.ylabel('data')
        plt.title(kwargs['title'])
        plt.show()


def create_quadratic_model(initial_values_dict):
    """
    Args:
        initial_values_dict: dict: {'a_0': value_1, 'a_1': value_2, 'a_2': value_3}

    Returns:
        params: instance of Parameters
        objective_function: objective function
    """
    params = Parameters()
    params.add('a_0', value=initial_values_dict['a_0'])
    params.add('a_1', value=initial_values_dict['a_1'], min=-10., max=10.)
    params.add('a_2', value=initial_values_dict['a_2'], min=-10., max=10.)

    def obj_func(params, x, y):
        # retrieve parameters of the model
        a_0 = params['a_0']
        a_1 = params['a_1']
        a_2 = params['a_2']

        # define exponentially decayed cosine with quadratic trend
        model = a_0 + a_1 * x + a_2 * x ** 2

        # return residuals
        return model - y

    return params, obj_func


def create_decayed_cosine_model():

    params = Parameters()
    params.add('amp', value=0.001, min=0)
    params.add('shift', value=0.0, min=-np.pi / 2., max=np.pi / 2)
    params.add('omega', value=0.3, min=0)
    params.add('decay', value=0.3, min=0)
    params.add('a_0', value=0., min=-10, max=10)
    params.add('a_1', value=0., min=-10, max=10)
    params.add('a_2', value=0., min=-10, max=10)

    # define objective function: it must return an array of residuals given params and the data x, y
    def obj_func(params, x, y):
        # retrieve parameters of the model
        amp = params['amp']
        shift = params['shift']
        omega = params['omega']
        decay = params['decay']
        a_0 = params['a_0']
        a_1 = params['a_1']
        a_2 = params['a_2']

        # define exponentially decayed cosine with quadratic trend
        model = a_0 + a_1 * x + a_2 * x ** 2 + amp * np.cos(x * omega + shift) * np.exp(- x * decay)

        # return residuals
        return model - y

    return params, obj_func


# Example of how to use the class NonLinearFitting:
def main_example():

    # create Parameters instance and add parameters to the model
    params = Parameters()
    params.add('a_0', value=0, min=0)
    params.add('a_1', value=0, min=0)
    params.add('a_2', value=0, min=0)
    params.add('decay', value=0.1)
    params.add('amp', value=0.1)
    params.add('shift', value=0.0, min=-np.pi / 2., max=np.pi / 2)
    params.add('omega', value=0.3, min=0)

    # define objective function: it must return an array of residuals given params and the data x, y
    def obj_func(params, x, y):
        # retrieve parameters of the model
        amp = params['amp']
        shift = params['shift']
        decay = params['decay']
        omega = params['omega']
        a_0 = params['a_0']
        a_1 = params['a_1']
        a_2 = params['a_2']

        # define exponentially decayed cosine with quadratic trend
        model = a_0 + a_1 * x + a_2 * x ** 2 + amp * np.cos(x * omega + shift) * np.exp(- x * decay)

        # return residuals
        return model - y

    # create data to be fitted
    x = np.linspace(0, 15, 301)
    y = 0.1 + 0.02 * x + 5. * np.cos(2 * x - 0.1) * np.exp(- x * 0.025) + np.random.normal(size=len(x), scale=0.2)

    # provide fitting
    non_linear_fitter = NonLinearFitting(obj_func, params, x, y)
    result = non_linear_fitter.least_squares()
    non_linear_fitter.plot_fit(title='Fitting exponentially decayed cosine with quadratic trend')

if __name__ == '__main__':
    pass
