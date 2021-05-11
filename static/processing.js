/*
# PBLMM - Peptide based linear mixed models for differential expression analysis
# Copyright (C) 2021 Kevin Klann
#This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

const {PythonShell} =require('python-shell')
const csv = require("csv-reader");
const fs = require("fs");
const {dialog} = require('electron').remote
const { ipcRenderer } = require('electron');
const { writableHighWaterMark } = require('csv-reader');
const Plotly = require("../static/plotly-latest.min.js");

function WaitForAnalysis(cb) {
    //Define button and interactive elements
    const out_plots=document.getElementById("output_plots");
    
    out_plots.style.display = "none"
    const start_button = document.getElementById("start");
    const file_button = document.getElementById("input_file");
    const design_button = document.getElementById("input_Design");
    file_button.addEventListener('click', function(){
        //On click open remote Renderer dialog box to get filepath for input    
        var file = dialog.showOpenDialog({ properties: ['openFile'] });
        file.then(function (result){
            window.file = result.filePaths[0];
            
        });
    });
    design_button.addEventListener('click', function(){
        //On click open remote Renderer dialog box to get filepath for design matrix    
        var design = dialog.showOpenDialog({ properties: ['openFile'] });
        design.then(function (result){
            window.design = result.filePaths[0];
        })
    });
    start_button.addEventListener("click", function(){
        if (window.result){
            clear()
            
        }
        //On click perform analysis
        console.log("Start")
        console.log(window.file)
        
        start_button.innerText = "Running..."
        //Get Parameters from text input fields
        const seq_col = document.getElementById("Sequence_column").value;
        const abun_col = document.getElementById("Abundance_column").value;
        const acc_col = document.getElementById("Accession_column").value;
        const norm = document.getElementById("Normalization").value;
        //create json with file paths and parameters
        let data1 = {
                "file1":window.file,
                "file2":window.design,
                "seq_col":seq_col,
                "abun_col":abun_col,
                "acc_col":acc_col,
                "norm":norm
                
            }
        //Invoke new python process with the Linear mixed models    
        let pyshell = new PythonShell('../app-1.0.4/resources/app/PBLMM.py',{pythonPath : 'python'});//../app-1.0.0/resources/app
        pyshell.send(JSON.stringify(data1),{mode:'json'});//Send data to python
        let messages = new Array();
        pyshell.on('message', function(message)  {
            messages.push(message); //Create array of all stdout stream messages
            console.log(JSON.parse(message)); //Returned stdout from python is stringifyed JSON
            
            window.result = JSON.parse(messages[messages.length - 1]); //messages.length should be 1, however if longer return only first element
            start_button.innerText = "Start analysis";//Reset button
            
        })
        pyshell.end(function (err) {
            //End python process
            if (err){
              throw err;
            };
            console.log("Finished");
            
            out_plots.style.display = "block";
            // Make plot element visible
            cb(); //callback to plot results


          });
    });
    };
function plot_p(){
    var plot_container_p = document.getElementById("pvalues");
        
        for (const condition of Object.keys(window.result)){
            let div = document.createElement('div');
            
            plot_container_p.appendChild(div);
            x1 = result[condition]['p_value'];
            var trace = {
                x: x1,
                type : 'histogram',
                name : 'P values',
                histnorm: "count",
               
                marker:{
                    color:"rgba(86, 211, 221, 0.5)",
                    
                    line:{
                        color:"rgba(86, 211, 221, 1)",
                        width:1
                    }
                },
                xbins:{
                    end: 1,
                    size: 0.01,
                    start: 0
                },
                autobinx:false
            }
            var data = [trace];
            var layout = {
                bargap: 0.5,
                barmode: "overlay",
                title: "P values "+String(condition),
                xaxis: {title:"P value"},
                yaxis: {title: "Count"}
            };
            Plotly.newPlot(div,data,layout) 
        };
}



function plot_fc(){
    var plot_container = document.getElementById("fcs");
        
        let div = document.createElement('div')
                     
        
        plot_container.appendChild(div)
        let traces = new Array()
        for (const condition of Object.keys(window.result)){
            
            x1 = result[condition]['fold_change']
            var trace = {
                y: x1,
                type : 'box',
                name : condition,
                points:'none',
                
                
                }
                
            traces.push(trace)    
            };
            
        
        var data =traces;
        var layout = {
            title:"Fold changes"
        }
        Plotly.newPlot(div,data,layout) 
};

function clear(){
    let pvalues = document.getElementById("pvalues");
    let fcs = document.getElementById("fcs");
    while (pvalues.lastElementChild) {
        pvalues.removeChild(pvalues.lastElementChild);
      }
    while (fcs.lastElementChild) {
        fcs.removeChild(fcs.lastElementChild);
    } 
    console.log("Cleared")
} 

function plot_results(){
    //WaitForAnalysis calls callback as soon as result is returned as JSON from the python script
    WaitForAnalysis(function(){
        plot_p();
        plot_fc();
    });
        //Get plot-container and append child row for each condition and plot Pvalues, qValues and fold changes
        
    
};



document.addEventListener("DOMContentLoaded", plot_results);//When site is rendered start script

//Install dependencies by python script
const deps = document.getElementById("dependencies");
deps.addEventListener("click", function (){
    let pyshell = new PythonShell('../app-1.0.4/resources/app/install_deps.py',{pythonPath : 'python3'})
    pyshell.on('message', function(printed)  {
        
        console.log(printed)
        
        
        
    })
    pyshell.end(function (err) {
        if (err){
          throw err;
        };
        options = {
            type : 'info',
            title : 'Install dependencies',
            message: "Finished"
        };
        dialog.showMessageBox(null,options)
        
      });
})
