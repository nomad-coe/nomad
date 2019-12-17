<#import "template.ftl" as layout>
<@layout.registrationLayout displayInfo=true; section>
    <#if section = "title">
        ${msg("updatePasswordTitle")?no_esc}
    <#elseif section = "header">
        ${msg("updatePasswordTitle")?no_esc}
    <#elseif section = "form">
        <form id="kc-passwd-update-form" class="form update-password ${properties.kcFormClass!}" action="${url.loginAction}" method="post">
            <input type="text" readonly value="this is not a login form" style="display: none;">
            <input type="password" readonly value="this is not a login form" style="display: none;">

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">lock</i>
                <input required id="password-new" class="mdc-text-field__input ${properties.kcInputClass!}" name="password-new" type="password">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="password" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("passwordNew")?no_esc}</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">lock</i>
                <input required id="password-confirm" class="mdc-text-field__input ${properties.kcInputClass!}" name="password-confirm" type="password">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="password-confirm" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("passwordConfirm")?no_esc}</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

            <div class="${properties.kcFormGroupClass!} update-password-button-container" >
                <div id="kc-form-buttons" class="${properties.kcFormButtonsClass!}">
                    <button class="mdc-button mdc-button--raised ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" type="submit">
                        ${msg("doSubmit")?no_esc}
                    </button>
                </div>
            </div>
        </form>
    </#if>
</@layout.registrationLayout>