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
const {app, BrowserWindow} = require('electron');
const { ipcMain } = require('electron');
const { ipcRenderer }=require('electron');
const dialog = require('electron').dialog;



function createWindow () {
    window = new BrowserWindow({width: 1800, height: 1000,frame: true, webPreferences: {
        enableRemoteModule: true,
        nodeIntegration: true,
        

    }})
    window.loadFile('./templates/main.html')
  }


app.on('ready', createWindow)

app.on('window-all-closed', function(){
    if(process.platform != 'darwin'){
        app.quit();
    }
});