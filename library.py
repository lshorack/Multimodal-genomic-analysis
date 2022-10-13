import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from sklearn.model_selection import train_test_split, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

font_title_defaults = {'fontsize': 25,
                   'fontweight': 'bold',
                   'color': 'black',
                   'alpha': 1,
                   'horizontalalignment': 'center',
                   'rotation': 0,
                   'family': 'georgia',
                   'style': 'italic'}

font_xlabel_defaults = {'fontsize': 12,
                   'fontweight': 'bold',
                   'color': 'black',
                   'alpha': 1,
                   'horizontalalignment': 'center',
                   'rotation': 0,
                   'family': 'georgia'}

font_ylabel_defaults = {'fontsize': 12,
                   'fontweight': 'bold',
                   'color': 'black',
                   'alpha': 1,
                   'horizontalalignment': 'center',
                   'rotation': 90,
                   'family': 'georgia'}


def fit_and_evaluate_citeseq_models(pipe, selected_genes, X_train, X_test, y_train, y_test, figures = False, pca_viz = True, eval_lr_coefs = False, cross_val = False):
    train_r_squared = []
    train_mean_squared_error = []
    train_mean_absolute_error = []
    
    test_r_squared = []
    test_mean_squared_error = []
    test_mean_absolute_error = []
    
    cross_validation_scores = {}
    lr_coefs = {}
    
    loops = 0
    for protein, gene_cols in zip(y_train.columns[4:], selected_genes): 
        # Fit the model.
        y_true = y_train[protein]      
        pipe.fit(X_train[gene_cols], y_true)
        y_preds = pipe.predict(X_train[gene_cols])
        
        # Evaluate Model with train data
        train_mean_squared_error.append(mean_squared_error(y_true, y_preds))
        train_r_squared.append(round(r2_score(y_true, y_preds), 2))
        train_mean_absolute_error.append(mean_absolute_error(y_true, y_preds))
        
        # Apply Cross Validation
        if cross_val:
            cross_validation_scores[protein] = cross_validate(pipe, X_train[selected_genes[0]], y_true, return_train_score = True)
        
        # Overwrite the true values with the test data.
        y_true = y_test[protein]
        y_preds = pipe.predict(X_test[gene_cols])
    
        # Evaluate Model with test data
        test_mean_squared_error.append(mean_squared_error(y_true, y_preds))
        test_r_squared.append(round(r2_score(y_true, y_preds), 2))
        test_mean_absolute_error.append(mean_absolute_error(y_true, y_preds))
        
        # Create Evaluatation Visualizations the model
        if figures:
            visualizations_linear_regression_model_evaluation(pipe, y_true, y_preds, protein, pca_viz)
        
        # Gather Coefs
        if eval_lr_coefs:
            lr_coefs[protein] = pd.DataFrame({'Coefficients': pipe['lr'].coef_, }, index = selected_genes[0])\
                                    .sort_values(by = 'Coefficients', ascending = False)

        if (loops % 14) == 0:
            print(f'{loops/140*100}% complete')
        loops += 1
    
    pipe_scores = pd.DataFrame({'Train R-Squared' : train_r_squared,
                              'Train Mean-Squared Error' : train_mean_squared_error,
                              'Train Mean-Absolute Error' : train_mean_absolute_error,
                              'Test R-Squared' : test_r_squared,
                              'Test Mean-Squared Error' : test_mean_squared_error,
                              'Test Mean-Absolute Error' : test_mean_absolute_error,
                             }, index = y_train.columns[4:])
    
    return pipe_scores, cross_validation_scores, lr_coefs


def visualizations_linear_regression_model_evaluation(pipe, y_true, y_preds, name, pca_viz = True):
    residuals = y_true - y_preds
    
    # Residual Plots
    fig, ax = plt.subplots(figsize = (12, 8))

    sns.residplot(x = y_preds,
                  y = residuals,
                  color = 'tab:cyan',
                  lowess = True, #help visualize relationship
                  line_kws = {'color': 'green'})

    ax.set_title('Residuals Plot', fontdict = font_title_defaults)
    ax.set_xlabel('Predicted Value', fontdict = font_xlabel_defaults)
    ax.set_ylabel('Residual Value', fontdict = font_ylabel_defaults)

    fig.savefig('./figures/residuals/' + 'residual_plot_' + name + '.png')
    plt.close()
    
    # PCA Plots
    if pca_viz:
        pca_var_ratio = pipe['pca'].explained_variance_ratio_

        # explained variance of components
        fig, ax = plt.subplots(figsize=(8,6))
        ax.plot(range(1,1001), pca_var_ratio, lw=2)
        ax.scatter(range(1, 1001), pca_var_ratio, s=100)
        ax.set_title('Explained Variance of Components', fontdict = font_title_defaults)
        ax.set_xlabel('Principal Component', fontdict = font_xlabel_defaults)
        ax.set_ylabel('Explained Variance', fontdict = font_ylabel_defaults);

        fig.savefig('./figures/pca_explained_variance/' + 'pca_exp_var_plot_' + name + '.png')
        plt.close()

        # cumulative explained variance of components
        cum_var_exp = np.cumsum(pca_var_ratio) * 100
        plt.figure(figsize=(9,7))
        component_number = range(1, 1001)
        plt.plot(component_number, cum_var_exp, lw=7)
        plt.axhline(y=0, linewidth=5, color='grey', ls='dashed')
        plt.axhline(y=100, linewidth=3, color='grey', ls='dashed')
        ax = plt.gca()
        ax.set_xlim([1,1001])
        ax.set_ylim([-5,105])
        ax.set_ylabel('Cumulative Variance Explained', fontdict = font_ylabel_defaults)
        ax.set_xlabel('Component', fontdict = font_xlabel_defaults)  
        ax.set_title('Component vs Cumulative Variance Explained\n', fontdict = font_title_defaults);

        fig.savefig('./figures/pca_cumulative_explained_variance/' + 'pca_cummulative_exp_var_plot_' + name + '.png')
        plt.close()