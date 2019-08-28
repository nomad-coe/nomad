<#import "template.ftl" as layout>
<@layout.registrationLayout displayInfo=true; section>
    <#if section = "title">
        ${msg("emailForgotTitle")?no_esc}
    <#elseif section = "header">
        ${msg("emailForgotTitle")?no_esc}
    <#elseif section = "form">
        <form id="kc-reset-password-form" class="form reset-password ${properties.kcFormClass!}" action="${url.loginAction}" method="post">
            <div class="reset-password-field ${properties.kcFormGroupClass!}">

                <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!} <#if usernameEditDisabled??>mdc-text-field--disabled</#if>">
                    <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">person</i>
                    <input required id="username" class="mdc-text-field__input ${properties.kcInputClass!}" name="username" type="text" autofocus>
                    <div class="${properties.kcLabelWrapperClass!}">
                        <label for="username" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("email")?no_esc}</label>
                    </div>
                    <div class="mdc-notched-outline">
                        <svg>
                            <path class="mdc-notched-outline__path"/>
                        </svg>
                    </div>
                    <div class="mdc-notched-outline__idle"></div>
                </div>

            </div>

            <div class="${properties.kcFormGroupClass!}">
                <div id="kc-form-buttons" class="${properties.kcFormButtonsClass!}">
                    <#--  <input class="btn btn-primary btn-flat btn-block ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" type="submit" value="${msg("doSubmit")}"/>  -->
                    <button class="mdc-button mdc-button--raised ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" type="submit">
                        ${msg("doSubmit")?no_esc}
                    </button>
                </div>
            </div>
            <div class="clearfix"></div>
        </form>
    <#elseif section = "info" >
        <hr />
        ${msg("emailInstruction")?no_esc}
    </#if>
</@layout.registrationLayout>
