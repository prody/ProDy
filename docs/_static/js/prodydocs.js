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

  // Alter logo based on URL
  var url = document.URL;
  var logo = document.getElementById("logo");
  if (url.search("plugins") > -1 || url.search("nmwiz") > -1) {
    logo.src = "http://www.csb.pitt.edu/ProDy/_static/nmwiz.png";
  } else if (url.search("evol") > -1 || url.search("panel1-3") > -1) {
    logo.src = "http://www.csb.pitt.edu/ProDy/_static/evol.png";
  }

  // Downloads
  $('tt.download').append('&nbsp;<i class="icon-download"></i>')
  var ext = $('a.external');
  ext.attr('target', '_blank');
  ext.each(function(i) {
    var jthis = $(this);
    var pre = jthis.find('span.pre');
    if (pre.length) {
      pre.append('&nbsp;<i class="icon-share"></i>');
    } else {
      jthis.append('&nbsp;<i class="icon-share"></i>');
    }
  });

});