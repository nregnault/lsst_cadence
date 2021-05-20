#!/usr/bin/env python 


"""
"""

import os
import os.path as op
import sys
import re
import glob

import pickle

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)

from argparse import ArgumentParser

import numpy as np

from yattag import Doc, indent

import cadence_io
from ubercal.cadence import Cadence



def dump_cadence_row(c, doc, root_dir='.'):
    """
    """
    tag, text = doc.tag, doc.text
    
    with tag('tr'):
        with tag('td'):
            text(c.name)
        with tag('td', style='text-align: center'):
            text('%10d' % c.nsn('faint'))
        with tag('td', style='text-align: center'):
            text('%5.2f' % c.zmax('faint'))
        with tag('td', style='text-align: center'):
            text('%10d' % c.nsn('normal'))
        with tag('td', style='text-align: center'):
            text('%5.2f' % c.zmax('normal'))
        with tag('td', style='text-align: center'):
            text('%5.3f' % c.cadence())            

def dump_cadence_summary_table(cads, doc, root_dir='.'):
    """
    """
    tag, text = doc.tag, doc.text
    
    with tag('table', id='myTable2', style='width:60%', align='center'):
        with tag('tr'):
            with tag('th', style='text-align: center', onclick="sortTable(0, 'myTable2')"):
                text('cadence name')
            with tag('th', style='text-align: center', onclick="sortTable(1, 'myTable2')"):
                text('NSN (faint)')
            with tag('th', style='text-align: center', onclick="sortTable(2, 'myTable2')"):
                text('zmax (faint)')
            with tag('th', style='text-align: center', onclick="sortTable(3, 'myTable2')"):
                text('nsn (normal)')
            with tag('th', style='text-align: center', onclick="sortTable(4, 'myTable2')"):
                text('zmax (normal)')
            with tag('th', style='text-align: center', onclick="sortTable(5, 'myTable2')"):
                text('median cadence')
                
        for c in cads.values():
            dump_cadence_row(c, doc, root_dir=root_dir)


def dump_cadence_summary_plots(cads, doc, root_dir='.'):
    """
    """
    tag, text = doc.tag, doc.text
    
    with tag('table', id='myTable3', style='width:100%', align='center'):
        with tag('tr'):
            with tag('th', style='text-align: center', onclick="sortTable(0, 'myTable3')"):
                text('cadence name')
            for k in ['seeing', 'survey_depth', 'compressed_cadence', 'compressed_nvisits', 'gap_10_compressed_cadence', 'gap_15_compressed_cadence']:                
                with tag('th', style='text-align: center'):
                    text(k)
                
        for c in cads.values():
            with tag('tr'):
                with tag('td'):
                    text(c.name)
                for k in ['seeing', 'survey_depth', 'compressed_cadence', 'compressed_nvisits', 'gap_10_compressed_cadence', 'gap_15_compressed_cadence']:
                    with tag('td'):
                        with tag('a', href=op.relpath(c.plots[k], root_dir)):
                            with tag('img', src=op.relpath(c.plots[k], root_dir), width='200px', style='text-align:center'):
                                pass
                

            
def get_html_head(doc, title=''):
    """
    """
    with doc.tag('head'):
        with doc.tag('title'):
            doc.text(title)
        doc.stag('meta', charset="UTF-8")
        doc.stag('meta', name="viewport", content="width=device-width, initial-scale=1")
        
        doc.stag('link', rel="stylesheet", href="https://www.w3schools.com/w3css/4/w3.css")
        doc.stag('link', rel="stylesheet", href="https://www.w3schools.com/lib/w3-theme-black.css")
        doc.stag('link', rel="stylesheet", href="https://fonts.googleapis.com/css?family=Roboto")
        doc.stag('link', rel="stylesheet", href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css")        
        doc.stag('link', rel="stylesheet", type="text/css", href="style.css")
        
        with doc.tag('style'):
            doc.text("""html,body,h1,h2,h3,h4,h5,h6 {font-family: "Roboto", sans-serif;}
            .w3-sidebar {
            z-index: 3;
            width: 250px;
            top: 43px;
            bottom: 0;
            height: inherit;
            }""")
                 
        with doc.tag('script', src='https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'):
            pass
        doc.asis('<script type="text/x-mathjax-config">MathJax.Hub.Config({TeX: { equationNumbers: { autoNumber: "AMS" } }});</script>')
        with doc.tag('script', src="http://html5shiv.googlecode.com/svn/trunk/html5.js"):
            pass
        with doc.tag('script', src='pager.js'):
            pass
        with doc.tag('script', src='sorttable.js'):
            pass

        
def get_navbar(doc, root_dir='.'):
    """
    """
    with doc.tag('div', klass='w3-top'):
        with doc.tag('div', klass='w3-bar w3-theme w3-top w3-left-align w3-large'):
            with doc.tag('a', klass='w3-bar-item w3-button w3-right w3-hide-large w3-hover-white w3-large w3-theme-l1', href="javascript:void(0)", onclick="w3_open()"):
                with doc.tag('i', klass='fa fa-bars'):
                    pass
            with doc.tag('a', href='index.html', klass='w3-bar-item w3-button w3-hide-small w3-hover-white'):
                doc.text('Summary')
            with doc.tag('a', href='pager.html', klass='w3-bar-item w3-button w3-hide-small w3-hover-white'):
                doc.text('Plots')
            with doc.tag('a', href='#', klass='w3-bar-item w3-button w3-hide-small w3-hover-white'):
                doc.text('About')
        

def get_sidebar(doc, cads, root_dir='.'):
    """
    """
    with doc.tag('nav', klass='w3-sidebar w3-bar-block w3-collapse w3-large w3-theme-l5 w3-animate-left', id='mySidebar'):
        with doc.tag('a', href='javascript:void(0)', onclick='w3_close()', klass='w3-right w3-xlarge w3-padding-large w3-hover-black w3-hide-large', title='Close Menu'):
            with doc.tag('i', klass='fa fa-remove'):
                pass
        with doc.tag('h4', klass='w3-bar-item'):
            with doc.tag('b'):
                doc.text('Menu')
        names = list(cads.keys())
        names.sort()
        for nm in names:
            with doc.tag('a', klass='w3-bar-item w3-button w3-hover-black', href=nm + '.html'):
                doc.text(nm)                
        with doc.tag('div', klass='w3-overlay w3-hide-large', onclick='w3_close()', style='cursor:pointer', title='close side menu', id='myOverlay'):
            pass


def get_pagination(doc, root_dir='.'):
    """
    """
    with doc.tag('div', klass='w3-center w3-padding-32'):
        with doc.tag('div', klass='w3-bar'):
            with doc.tag('a', klass='w3-button w3-black', href='#'):
                doc.text(1)
            

def get_footer(doc, root_dir='.'):
    """
    """
    with doc.tag('footer', id='myFooter'):
        with doc.tag('div', klass='w3-container w3-theme-l2 w3-padding-32'):
            with doc.tag('h4'):
                doc.text('Footer')
        with doc.tag('div', klass='w3-container w3-theme-l1'):
            with doc.tag('p'):
                doc.text('Powered by ')
                with doc.tag('a', href='https://www.w3schools.com/w3css/default.asp', target='_blank'):
                    doc.text('w3.css')
                                    
        
def html_post(cads, root_dir='.'):
    """
    """
    doc, tag, text = Doc().tagtext()

    doc.asis('<!DOCTYPE html>')
    
    with tag('html'):
        get_html_head(doc, title='Analysis of the fbs 4.1 cadences')

        
        with tag('body'):
            get_navbar(doc, root_dir=root_dir)
            get_sidebar(doc, cads, root_dir=root_dir)
            with tag('div', klass='w3-main', style='margin-left:250px'):
                with tag('div', klass='w3-row w3-padding-64'):
                    with tag('div', klass='w3-container'):

                        with tag('h1', style='text-align:center'):
                            text('Introduction')

                        with tag('h1', style='text-align:center'):
                            text('Summary metrics')
                            
                        dump_cadence_summary_table(cads, doc, root_dir=root_dir)

                        with tag('h1', style='text-align:center'):
                            text('Global plots')
                        dump_cadence_summary_plots(cads, doc, root_dir=root_dir)
                        
            get_pagination(doc, root_dir=root_dir)
            get_footer(doc, root_dir=root_dir)
            
    return indent(doc.getvalue())



def pager_page(cads, root_dir='.', plot_list=['seeing.png', 'survey_depth.png',
			                      'compressed_cadence.png', 'compressed_nvisits.png',
			                      'raw_nvisits.png', 'gap_10_compressed_cadence.png', 'gap_15_compressed_cadence.png', ]):
    """
    """
    doc, tag, text = Doc().tagtext()

    doc.asis('<!DOCTYPE html>')


    with tag('html'):
        get_html_head(doc, title='Pager')

        with tag('body'):
            get_navbar(doc, root_dir=root_dir)
            sorted_cadence_names = [c.name for c in cads.values()]
            sorted_cadence_names.sort()
            with tag('script'):
                text("""
                var cadence = {};
                cadence.name = 'cadence';
                cadence.value = 'none';
                var sn = {};
                sn.name = 'sn';
                sn.value = 'normal';
                var plot = {};
                plot.name = 'plot';
                plot.value = 'seeing.png';
                var cadences = [""" + ', '.join(['"' + nm + '"' for nm in sorted_cadence_names]) + '];')

            # adapted sidebar
            with tag('nav', klass='w3-sidebar w3-bar-block w3-collapse w3-large w3-theme-l5 w3-animate-left', id='mySidebar'):
                with tag('a', href='javasrcript:void(0)', onclick='w3_close()', klass='w3-right w3-xlarge w3-padding-large w3-hover-black w3-hide-large', title='Close Menu'):
                    with tag('i', klass='fa fa-remove'):
                        pass
                with tag('script'):
                    text('make_menu(cadences, cadence, 100, "left");')
                    
            # main body
            with tag('div', klass='w3-main', style='margin-left:250px'):
                with tag('div', klass='w3-row w3-padding-64'):
                    with tag('div', klass='w3-container'):
                        
                        with tag('div', id='left', style='width:400px; float:left;'):
                            with tag('div', id='menus', style='width:100%'):
                                with tag('script'):
                                    text("""var sne= ["normal", "faint"];
                                    make_menu(sne, sn, 100, 'left');
                                    var plots = %r;
                                    make_menu(plots, plot, 200, 'left');""" % plot_list)
                            with tag('table', id='more_buttons', style='width:100%; margin-left:50px; margin-top:50px; display:inline-block;'):
                                with tag('tr'):
                                    with tag('td', klass='menu_entry'):
                                        with tag('a', klass='menu_anchor', href='javascript:next_cadence_same_plot();'):
                                            text('Next cadence same plot')
                                with tag('tr'):
                                    with tag('td', klass='menu_entry'):
                                        with tag('a', klass='menu_anchor', href='javascript:next_plot_same_cadence();'):
                                            text('Next plot same cadence')
                        with tag('div', klass='w3-center w3-large'):
                            with tag('h1', id='figure_cadence_name'):
                                pass
                            with tag('img', id='figure', src='none', widht='1000px', style='float:center'):
                                pass
                    
                with tag('script'):
                    text("""set_value(cadence, cadences[0]);
                            set_value(sn, sne[0]); 
                            set_value(plot, plots[0]); 
                            assign_plot();""")
                             
        get_pagination(doc, root_dir=root_dir)
        get_footer(doc, root_dir=root_dir)
                    
    return indent(doc.getvalue())


def cadence_summary_page(c, cads, root_dir='.'):
    """
    """

    doc, tag, text = Doc().tagtext()

    doc.asis('<!DOCTYPE html>')    
    with tag('html'):
        get_html_head(doc, title='Metrics for cadence: %s' % c.name)
        
        with tag('body'):
            get_navbar(doc, root_dir=root_dir)
            get_sidebar(doc, cads, root_dir=root_dir)
            with tag('div', klass='w3-main', style='margin-left:250px'):
                with tag('div', klass='w3-row w3-padding-64'):
                    with tag('div', klass='w3-container'):

                        with tag('div', style='text-align:center'):
                            with tag('h1'):
                                text('Cadence: %s' % c.name)

                        with tag('h2', style='text-align:center', background='yellow'):
                            text('Normal SN')
                            
                        with tag('div', style='text-align:center'):
                            with tag('video', width='80%', controls=1):
                                with tag('source', src=op.relpath(c.movie_file['normal'], root_dir)):
                                    text('Your browser does not support HTML5 videos')

                        with tag('h2', style='text-align:center', background='yellow'):
                            text('Faint SN')
                            
                        with tag('div', style='text-align:center'):
                            with tag('video', width='80%', controls=1):
                                with tag('source', src=op.relpath(c.movie_file['faint'], root_dir)):
                                    text('Your browser does not support HTML5 videos')

                                    
            get_pagination(doc, root_dir=root_dir)
            get_footer(doc, root_dir=root_dir)
            
    return indent(doc.getvalue())


def generate_report(cads, root_dir='html_report'):
    """
    """
    if not op.isdir(root_dir):
        os.makedirs(root_dir)

    #    old_cwd = os.getcwd()

    try:
        #       os.chdir(root_dir)
        
        # dump the cadence data 
        cc = cadence_io.dump_cadences(cads, root_dir + os.sep + 'plots')
        cads = cc
    
        #    cads = load_cadences()
        p = html_post(cc, root_dir=root_dir)

        with open(root_dir + os.sep + 'index.html', 'w') as f:
            f.write(p)

        for c in cads.values():
            p = cadence_summary_page(c, cads, root_dir=root_dir)
            with open(root_dir + os.sep + c.name + '.html', 'w') as f:
                f.write(p)
        # pager
        pager = pager_page(cads)
        with open(root_dir + os.sep + 'pager.html', 'w') as f:
            f.write(pager)
                
    finally:
        # os.chdir(old_cwd)
        pass


if __name__ == "__main__":
    """
    """
    parser = ArgumentParser(description='')
    parser.add_argument('-O', '--output-dir', 
                        default='./', 
                        help='output directory')
    parser.add_argument('pickled_db',
                        help='input pickled cadence database')

    args = parser.parse_args()

    #    with open(sys.argv[1], 'rb') as f:
    #        cads = pickle.load(f)
    with open(args.pickled_db, 'rb') as f:
        cads = pickle.load(f)

    generate_report(cads, root_dir=args.output_dir)


    
