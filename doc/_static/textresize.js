$(document).ready(function() {

  fsUnit = 'px';
  defaultFS = '12';
  function changeFontSize(incr) {
    if (document.body.style.fontSize == "") {
      document.body.style.fontSize = defaultFS + fsUnit;
    }
    if (incr == 0) {
      document.body.style.fontSize = defaultFS + fsUnit;
    } else {
      document.body.style.fontSize = (parseFloat(document.body.style.fontSize) + incr ) + fsUnit;
    }
    var exdate = new Date();
    exdate.setDate(exdate.getDate() + 365);
    document.cookie= "ProDyDocFontSize" + "=" + escape(parseFloat(document.body.style.fontSize)) + ";expires=" + exdate.toUTCString();
  }

  function setFontSize(fs) {
    document.body.style.fontSize = fs + fsu;
    var exdate = new Date();
    exdate.setDate(exdate.getDate() + 365);
    document.cookie= "ProDyDocFontSize" + "=" + escape(fs) + ";expires=" + exdate.toUTCString();
  }

  function getFontSize() {
    if (document.cookie.length > 0) {
      c_name = "ProDyDocFontSize"
      var c_start = document.cookie.indexOf(c_name + "=");
      if (c_start != -1) {
        c_start = c_start + c_name.length + 1;
        var c_end = document.cookie.indexOf(";", c_start);
        if (c_end == -1) c_end = document.cookie.length;
          return unescape(document.cookie.substring(c_start, c_end));
      }
    }
    return "";
  }

  function checkFontSize() {
    var fs = getFontSize('SphinxFontSize');
    if (fs != null && fs != "") {
      document.body.style.fontSize = fs + fsUnit;
    }
  }
  $('#text-larger').click(function() {changeFontSize(1)});
  $('#text-smaller').click(function() {changeFontSize(-1)});
  $('#text-reset').click(function() {changeFontSize(0)});
  checkFontSize();

});