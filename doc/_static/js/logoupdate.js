$(document).ready(function() {

    function updateLogo() {
      if (document.URL.search("\/plugins\/") > -1 || document.URL.search("panel1-4") > -1) {
        document.getElementById("logo").src = "http://www.csb.pitt.edu/ProDy/_static/nmwiz.png";
      } else if (document.URL.search("evol") > -1 || document.URL.search("panel1-3") > -1) {
        document.getElementById("logo").src = "http://www.csb.pitt.edu/ProDy/_static/evol.png";
      }

    }
    updateLogo();
});
