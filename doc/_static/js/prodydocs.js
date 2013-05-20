
$(document).ready(function() {

  var $window = $(window)
  var items = $("#localtoc a")
  if (items.length) {
    items.prepend('<i class="icon-chevron-right pull-right"></i>')
    items[0].href = "#" + $(".headerlink")[0].href.split("#")[1]
    $("#localtoc ul").addClass("nav nav-list")
    $("#localtoc > ul").addClass("bs-docs-sidenav")
  }


  setTimeout(function () {
      $('.bs-docs-sidenav,.subnav').affix({
        offset: {
          top: 140,
          bottom: 20
        }
      })
    }, 100);

});