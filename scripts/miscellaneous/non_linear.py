
from lmfit import Minimizer, Parameters


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
