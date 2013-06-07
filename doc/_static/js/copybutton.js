
$(document).ready(function() {
    /* Add a [>>>] button on the top-right corner of code samples to hide
     * the >>> and ... prompts and the output and thus make the code
     * copyable. */
    var div = $('.highlight-python .highlight,' +
                '.highlight-python3 .highlight,' +
                '.highlight-ipython .highlight')
    var pre = div.find('pre');
    // get the styles from the current theme
    pre.parent().parent().css('position', 'relative');
    var hide_text = 'Hide the prompts and output';
    var show_text = 'Show the prompts and output';
    var border_width = pre.css('border-top-width');
    var border_style = pre.css('border-top-style');
    var border_color = pre.css('border-top-color');
    var button_styles = {
        'cursor':'pointer', 'position': 'absolute', 'top': '0', 'right': '0',
        'border-color': border_color, 'border-style': border_style,
        'border-width': border_width, 'color': border_color, 'text-size': '75%',
        'font-family': 'monospace', 'padding-left': '0.2em', 'padding-right': '0.2em',
        'border-radius': '0 3px 0 0'
    }

    if (!(div.length && pre.find('.gp'))) {
        $("#showcodebuttons").hide();
    } else {
        var prepcode = function(wc) {
            var codesnippets = $('#codesnippets')
            codesnippets.html('')
            //codesnippets.append("#!/usr/bin/env python\n")
            //codesnippets.append("# -*- coding: utf-8 -*-\n")
            codesnippets.append("# This code was copied from ProDy documentation.\n")
            codesnippets.append("# Title: " + document.title.replace('\xe2', '-') + "\n")
            codesnippets.append("# URL: " + document.URL + "\n\n")
            pre.each(function(index) {
                jthis = $(this).clone();
                if (wc) {
                    jthis.find('.go:empty').remove();
                    jthis.find('.go').prepend('# ');
                } else {
                    jthis.find('.go').remove();
                }
                jthis.find('.gr, .gp').remove();
                var lines = jthis.text().split('\n');
                $.each(lines, function(l) {
                    var line = lines[l];
                    if (line.length) {
                        codesnippets.append(line.replace(/\</g,"&lt;").replace(/\>/g,"&gt;").replace(/\&/g,"&amp;").replace(/\'/g,"&apos;").replace(/\"/g,"&quot;") + '\n');
                    }
                });
                codesnippets.append('\n');
            });
        }

        $("#showcode").click(function () {
            prepcode(0);
        });
        $("#withcomments").click(function () {
            prepcode(1);
        });
    }

    $("#selectcode").click(function() {
        var e = document.getElementById('codesnippets')
        e.focus();
        e.select();
    })

    var clip = new ZeroClipboard($("#copycode"), {moviePath: $("#zeroclipboardpath").val()});

    // create and add the button to all the code blocks that contain >>>
    div.each(function(index) {
        var jthis = $(this);
        var gp = jthis.find('.gp');
        if (gp.length > 0) {
            if (jthis.parent().hasClass('highlight-ipython')) {
                var button = $('<span class="copybutton">' + gp[0].innerHTML.split(':')[0] + '</span>');
            }
            else {
                var button = $('<span class="copybutton">&gt;&gt;&gt;</span>');
            }
            button.css(button_styles)
            button.attr('title', hide_text);
            button.data("clicks", 1)
            jthis.prepend(button);
        }
        // tracebacks (.gt) contain bare text elements that need to be
        // wrapped in a span to work with .nextUntil() (see later)
        jthis.find('pre:has(.gt)').contents().filter(function() {
            return ((this.nodeType == 3) && (this.data.trim().length > 0));
        }).wrap('<span>');
    });

    // define the behavior of the button when it's clicked
    $('.copybutton').click(
        function() {
            var button = $(this);
            clicks = button.data('clicks');
            if (clicks) {
                button.parent().find('.go, .gp, .gt, .gr').hide();
                button.next('pre').find('.gt, .gr').nextUntil('.gp, .go').css('visibility', 'hidden');
                button.css('text-decoration', 'line-through');
                button.attr('title', show_text);
            } else {
                button.parent().find('.go, .gp, .gt, .gr').show();
                button.next('pre').find('.gt, .gr').nextUntil('.gp, .go').css('visibility', 'visible');
                button.css('text-decoration', 'none');
                button.attr('title', hide_text);
            }
            $(this).data("clicks", !clicks);
        });
});

