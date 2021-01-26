const {app, BrowserWindow} = require('electron');
const { ipcMain } = require('electron');
const { ipcRenderer }=require('electron');
const dialog = require('electron').dialog;



function createWindow () {
    window = new BrowserWindow({width: 1800, height: 1000,frame: false, webPreferences: {
        enableRemoteModule: true,
        nodeIntegration: true,
        devTools: false

    }})
    window.loadFile('./templates/main.html')
  }


app.on('ready', createWindow)

app.on('window-all-closed', function(){
    if(process.platform != 'darwin'){
        app.quit();
    }
});