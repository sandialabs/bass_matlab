[![Pipeline Status](https://github.com/sandialabs/BASS_matlab/actions/workflows/matlab.yml/badge.svg)](https://github.com/sandialabs/BASS_matlab/actions/workflows/matlab.yml)

# BASS Matlab

Implementation of Bayesian Adaptive Spline Surfaces in MATLAB

# References

Francom, Devin, Bruno Sansó, Vera Bulaevskaya, Donald Lucas, and Matthew
Simpson. 2019. “Inferring Atmospheric Release Characteristics in a Large
Computer Experiment Using Bayesian Adaptive Splines.” *Journal of the
American Statistical Association*.
<https://doi.org/10.1080/01621459.2018.1562933>.

Francom, Devin, Bruno Sansó, Ana Kupresanin, and Gardar Johannesson.
2018. “<span class="nocase">Sensitivity analysis and emulation for
functional data using Bayesian adaptive splines</span>.” *Statistica
Sinica*. <https://doi.org/10.5705/ss.202016.0130>.

## Using BASS inside a larger Gibbs sampler

`bass` fits a fixed response from start to finish. When the response is itself a latent variable updated by another part of a sampler (data augmentation), use `bassGibbs`. It keeps the chain's state between visits and advances it a few steps at a time against the current response:

```matlab
g = bassGibbs(X, y0, nstore, 'g1', 2, 'g2', 0.01);
for it = 1:n_iter
    % ... update the latent response y ...
    g.setResponse(y);   % keep the basis functions, swap the response
    g.step(3);          % 3 RJMCMC moves, each followed by Gibbs draws of s2, beta, lam, tau
    mu = g.fitted();    % regression function at the training inputs
    s2 = g.s2();
    if keep, g.save(); end
end
mod = g.model();        % BassModel with the saved draws: predict(), plot()
```

Options added for this use, available in both `bass` and `bassGibbs`. The defaults reproduce the original model exactly:

| option | default | meaning |
|---|---|---|
| `intercept` | `true` | include a constant basis function |
| `center_basis` | `false` | center each basis function over the training inputs. With `intercept = false`, every posterior draw of the regression function averages to exactly zero over the training inputs. |
| `s2_fixed` | `NaN` | hold the error variance at this value instead of sampling it |
