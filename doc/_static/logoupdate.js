$(document).ready(function() {

    function updateLogo() {
      if (document.URL.search("\/plugins\/") > -1 || document.URL.search("panel1-4") > -1) {
        document.getElementById("logo").src = "/_static/nmwiz.png";
      } else if (document.URL.search("evol") > -1 || document.URL.search("panel1-3") > -1) {
        document.getElementById("logo").src = "/_static/evol.png";
      }

    }
    updateLogo();
});