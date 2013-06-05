
$(document).ready(function() {

  var $window = $(window)
  var items = $("#localtoc a")
  if (items.length) {
    items.prepend('<i class="icon-chevron-right pull-right"></i>')

    var sections = $(".section");
    var section_id = sections[0].id;
    if (section_id.indexOf('.') >= 0) {
      items[0].href = "#" + $(sections[0]).find('>:first-child').attr('id');
    } else {
      items[0].href = "#" + section_id;
    }

    items.first().addClass("navfirst")
    items.last().addClass("navlast")
    $("#localtoc ul").addClass("nav nav-list")
    $("#localtoc > ul").addClass("toc")
  }

  setTimeout(function() {
      $('#sidebar,.subnav').affix({
        offset: {
          top: 140,
          bottom: 20
        }
      });
    }, 100);

  $('.collapse').not('#collapseToolbox').collapse('hide');

  $('#carousel').carousel({interval: 100000}).on({
      slid: function(){
        var pills = $('.nav-pills > li');
        pills.removeClass('active');
        var index = $('#carousel .item.active').index('#carousel .item');
        $(pills[index]).addClass('active');
        if (index == 1) {
          document.getElementById("logo").src = "_static/evol.png";
        } else if (index == 2) {
          document.getElementById("logo").src = "_static/nmwiz.png";
        } else {
          document.getElementById("logo").src = "_static/logo.png";
        }

      }});

  if (document.URL.search("plugins") > -1 || document.URL.search("nmwiz") > -1) {
    document.getElementById("logo").src = "http://www.csb.pitt.edu/ProDy/_static/nmwiz.png";
  } else if (document.URL.search("evol") > -1 || document.URL.search("panel1-3") > -1) {
    document.getElementById("logo").src = "http://www.csb.pitt.edu/ProDy/_static/evol.png";
  }

});