window.onload = function() {
    // add ripple effect to all material buttons on the page
    document.querySelectorAll('.mdc-button').forEach(function(e) {
        mdc.ripple.MDCRipple.attachTo(e);
    });

    // initialize all text fields
    document.querySelectorAll('.mdc-text-field').forEach(function(e) {
        new mdc.textField.MDCTextField(e);
    });

    // initialize all icons
    document.querySelectorAll('.mdc-text-field__icon').forEach(function(e) {
        new mdc.textField.MDCTextFieldIcon(e);
    });

    // initialize the language select box
    try {
        var select = new mdc.select.MDCSelect(
            document.querySelector('.language-picker .mdc-select')
        );

        select.listen('change', function() {
            var redirectUrl = document.querySelector('#language-picker-dropdown')
                .value;
            window.location.href = redirectUrl;
        });
    } catch {}

    const queryString = window.location.search
    const urlParams = new URLSearchParams(queryString)
    const redirectUri = urlParams.get('redirect_uri')
    const backLink = this.document.querySelector('#mdc-back-link')
    backLink.onclick = function() {
        if (redirectUri) {
            window.location.href = redirectUri;
        } else {
            window.history.back()
        }
    }
    backLink.innerHTML = redirectUri ? 'back to NOMAD' : 'back';
};
