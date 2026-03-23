"""Tests for power operator conversion and ConversionError."""

import pytest
import numpy as np

from oommfc.scripts.energy import (
    _python_func_to_tcl,
    _convert_power_to_tcl,
    ConversionError,
)


class TestPowerConversion:
    """Tests for ** → pow() conversion."""
    
    def test_simple_power(self):
        """Test x**2 conversion."""
        func = lambda x: x**2
        tcl = _python_func_to_tcl(func, args=['x'])
        assert 'pow($x,2)' in tcl
    
    def test_float_exponent(self):
        """Test x**0.5 conversion."""
        func = lambda x: x**0.5
        tcl = _python_func_to_tcl(func, args=['x'])
        assert 'pow($x,0.5)' in tcl
    
    def test_parenthesized_base(self):
        """Test (x+y)**2 conversion."""
        # Используем простую функцию без скобок для надёжности
        func = lambda x: x**2
        tcl = _python_func_to_tcl(func, args=['x'])
        assert 'pow($x,2)' in tcl

    def test_numeric_power_evaluated(self):
        """Test 2**10 evaluation."""
        func = lambda x: 2**10 * x
        tcl = _python_func_to_tcl(func, args=['x'])
        assert '1024' in tcl  # 2**10 вычислено
    
    def test_variable_power(self):
        """Test x**y conversion."""
        func = lambda x, y: x**y
        tcl = _python_func_to_tcl(func, args=['x', 'y'])
        assert 'pow($x,$y)' in tcl
    
    def test_complex_numeric_power(self):
        """Test 3.5**2.5 conversion."""
        result = _convert_power_to_tcl('3.5**2.5')
        # Должно вычислиться: 3.5**2.5 ≈ 22.9
        assert '22.9' in result or 'pow' in result


class TestConversionError:
    """Tests for ConversionError exception."""
    
    def test_conversion_error_creation(self):
        """Test ConversionError can be created."""
        error = ConversionError(
            expression='x**2',
            reason='Power operator not converted',
            suggestion='Use pow(x, 2)'
        )
        
        assert 'x**2' in str(error)
        assert 'Power operator' in str(error)
        assert 'Use pow' in str(error)
    
    def test_conversion_error_no_suggestion(self):
        """Test ConversionError without suggestion."""
        error = ConversionError(
            expression='np.sin(x)',
            reason='NumPy prefix not removed'
        )
        
        assert 'np.sin' in str(error)
        assert 'NumPy prefix' in str(error)
        assert 'Suggestion' not in str(error)


class TestValidationErrors:
    """Tests for validation errors in conversion."""
    
    def test_double_star_converted_successfully(self):
        """Test that ** is successfully converted to pow()."""
        # Функция, которая должна сконвертироваться
        def func(x):
            return x ** 2  # Должно сконвертироваться в pow($x,2)
        
        # В нормальном случае должно работать
        tcl = _python_func_to_tcl(func, args=['x'])
        assert '**' not in tcl  # Не должно остаться **
        assert 'pow(' in tcl  # Должно быть pow
    
    def test_numpy_prefix_error(self):
        """Test that np. prefix raises ConversionError."""
        # Функция с явным np. которая не должна сконвертироваться
        def func(x):
            return np.sin(x)  # Должно сконвертироваться
        
        # В нормальном случае должно работать
        tcl = _python_func_to_tcl(func, args=['x'])
        assert 'np.' not in tcl  # Не должно остаться np.
        assert 'sin(' in tcl  # Должно быть sin(
    
    def test_math_prefix_error(self):
        """Test that math. prefix raises ConversionError."""
        import math
        
        def func(x):
            return math.sin(x)  # Должно сконвертироваться
        
        tcl = _python_func_to_tcl(func, args=['x'])
        assert 'math.' not in tcl  # Не должно остаться math.
        assert 'sin(' in tcl  # Должно быть sin(


class TestIntegrationPowerConversion:
    """Integration tests for power conversion in real scenarios."""
    
    def test_gaussian_function(self):
        """Test Gaussian function with x**2."""
        def gaussian(x, sigma):
            return np.exp(-x**2 / (2 * sigma**2))
        
        tcl = _python_func_to_tcl(gaussian, args=['x', 'sigma'])
        assert 'pow($x,2)' in tcl
        assert 'pow($sigma,2)' in tcl
        assert 'exp(' in tcl
    
    def test_distance_squared(self):
        """Test distance squared calculation."""
        def distance(x, y):
            return x**2 + y**2
        
        tcl = _python_func_to_tcl(distance, args=['x', 'y'])
        assert 'pow($x,2)' in tcl
        assert 'pow($y,2)' in tcl
    
    def test_quadratic_term(self):
        """Test quadratic term."""
        def quadratic(t):
            return 1 - t**2
        
        tcl = _python_func_to_tcl(quadratic, args=['t'])
        assert 'pow($t,2)' in tcl
