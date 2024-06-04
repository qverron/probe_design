import os
import polars as pl
import lets_plot
from lets_plot import *
LetsPlot.setup_html()

# tested with 
# lets_plot==4.3.2
# polars==0.20.30

def read_summary(fname:str='data/final_probes_summary.tsv')-> pl.DataFrame:
    """
    Read the summary of the probes tsv file into a polars DataFrame
    default file is 'data/final_probes_summary.tsv'
    converts the 'roi' column to string

    Parameters
    ----------
    fname : str, optional
        file name of the summary, by default 'data/final_probes_summary.tsv'

    Returns
    -------
    pl.DataFrame
        DataFrame with the summary of the probes

    """

    summary = pl.read_csv(fname,
                        has_header=True,
                        separator='\t').with_columns(pl.col('roi').cast(pl.String))

    return summary

def add_stds(summary: pl.DataFrame)-> pl.DataFrame:
    """
    add (and subtract) the standard deviations to the summary DataFrame
    creating min and max values for the error bars

    Parameters
    ----------
    summary : pl.DataFrame
        DataFrame with the summary of the probes

    Returns
    -------
    pl.DataFrame
    """

    summary = summary.with_columns(
        d_mean_minus_std = pl.col('d_mean') - pl.col('d_std'),
        d_mean_plus_std = pl.col('d_mean') + pl.col('d_std'),
        gc_mean_minus_std = pl.col('gc_mean') - pl.col('gc_std'),
        gc_mean_plus_std = pl.col('gc_mean') + pl.col('gc_std'),
        mean_oligo_cost_minus_std = pl.col('mean_oligo_cost') - pl.col('std_oligo_cost'),
        mean_oligo_cost_plus_std = pl.col('mean_oligo_cost') + pl.col('std_oligo_cost'),
    )

    return summary


def plot_coverage(summary:pl.DataFrame)-> lets_plot.plot.core.PlotSpec:
    """
    Plot the coverage by ROI
    default palette is hue
    default columns are x='roi', y='coverage'

    summary: pl.DataFrame
        DataFrame with the summary of the probes

    Returns
    -------
    lets_plot.plot.core.PlotSpec
        Plot of the coverage by ROI

    """
    coverage_plot = (
        ggplot(summary) +
        geom_bar(aes(x='roi', y='coverage',fill='roi'),
                color='#3F3F3F', 
                stat='identity',
                tooltips=layer_tooltips().title('ROI: @roi').
                                            line('@chr').
                                            line('Coverage: @coverage').
                                            line('Probe Size: @probe_size').
                                            line('Region Size: @region_size')
                )+
        ggtitle('Coverage by Rois')+
        scale_fill_hue()+
        theme_minimal2()+
        theme(line=element_line(size=2,color='black'),
            #panel_grid=element_line(size=1,color='#AFAFAF'),
            legend_position='right',
            )+
        xlab('ROI')+
        ylab('Coverage')
    )+ggsize(800, 600)
    coverage_plot

    return coverage_plot


def plot_distance(summary: pl.DataFrame)-> lets_plot.plot.core.PlotSpec:
    """Draw a plot of the distance between oligos by ROI
    default palette is inferno

    Parameters
    ----------
    summary : pl.DataFrame
       summary of the probes

    Returns
    -------
    lets_plot.plot.core.PlotSpec
        Plot of the distance between oligos by ROI
    """


    distance_plot = (
        ggplot(summary) +
        geom_bar(aes(x='roi', y='d_mean',fill='roi'),
                color='#3F3F3F', 
                stat='identity',
                tooltips=layer_tooltips().title('ROI: @roi').
                                            line('@chr').
                                            line('Distance Mean: @d_mean').
                                            line('Probe Size: @probe_size')
                )+
        geom_errorbar(
                        aes(x='roi', ymin='d_mean_minus_std', ymax='d_mean_plus_std'),color='#3f3f3f',
                        tooltips=layer_tooltips().line('Std: @d_std')
                                                

                    )+
        ggtitle('Distance between oligos')+
        scale_fill_viridis(option='inferno')+
        theme_minimal2()+
        theme(line=element_line(size=2,color='black'),
            #panel_grid=element_line(size=1,color='#AFAFAF'),
            legend_position='right',
            )+
        xlab('ROI')+
        ylab('Distance')
        

    )+ggsize(800, 600)

    return distance_plot


def plot_oligo_cost(summary: pl.DataFrame)-> lets_plot.plot.core.PlotSpec:
    """
    Plot oligo cost by ROI
    default palette is viridis
    default columns are x='roi', y='mean_oligo_cost'
    
    summary: pl.DataFrame
        DataFrame with the summary of the probes\
    
    Returns
    -------
    lets_plot.plot.core.PlotSpec
        Plot of oligo cost by ROI
    
    """
    oligo_cost_plot = (
        ggplot(summary) +
        geom_bar(aes(x='roi', y='mean_oligo_cost',fill='roi'),
                color='#3F3F3F', 
                stat='identity',
                tooltips=layer_tooltips().title('ROI: @roi').
                                            line('@chr').
                                            line('Mean oligo cost: @mean_oligo_cost').
                                            line('Probe Size: @probe_size')
                )+
        geom_errorbar(
                        aes(x='roi', ymin='mean_oligo_cost_minus_std', ymax='mean_oligo_cost_plus_std'),color='#3f3f3f',
                        tooltips=layer_tooltips().line('Std: @std_oligo_cost') 
                    )+
        ggtitle('Oligo Costs')+
        scale_fill_viridis()+
        theme_minimal2()+
        theme(line=element_line(size=2,color='black'),
            legend_position='right',
            )+
        xlab('ROI')+
        ylab('Oligo Cost')
    )+ggsize(800, 600)

    return oligo_cost_plot


def main(filename, save_dir):
    summary = read_summary(filename) # read the summary file
    summary = add_stds(summary) # add the standard deviations to the summary
    coverage_plot = plot_coverage(summary) # plot the coverage
    distance_plot = plot_distance(summary) # plot the distance 
    oligo_cost_plot = plot_oligo_cost(summary) # plot the oligo cost

    # parent directory should exist before creating the save directory
    if not os.path.exists(save_dir):
        os.mkdir(save_dir) # create the save directory

    # html (interactive plots)
    coverage_plot.to_html(os.path.join(save_dir,'coverage_plot.html')) # save the coverage plot
    distance_plot.to_html(os.path.join(save_dir,'distance_plot.html')) # save the distance plot
    oligo_cost_plot.to_html(os.path.join(save_dir,'oligo_cost_plot.html')) # save the oligo cost plot
    

    # svg (vector graphics)
    coverage_plot.to_svg(os.path.join(save_dir,'coverage_plot.svg'))
    distance_plot.to_svg(os.path.join(save_dir,'distance_plot.svg'))
    oligo_cost_plot.to_svg(os.path.join(save_dir,'oligo_cost_plot.svg'))


    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Visual Report')
    parser.add_argument('-f','--filename', type=str, default='data/final_probes_summary.tsv', help='Summary file')
    parser.add_argument('-s','--save_dir', type=str, default='data/visual_summary', help='save (output) folder')
    args = parser.parse_args()
    main(filename=args.filename,save_dir=args.save_dir)