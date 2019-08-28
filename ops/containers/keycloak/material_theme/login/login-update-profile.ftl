<#import "template.ftl" as layout>
<@layout.registrationLayout; section>
    <#if section = "title">
        ${msg("loginProfileTitle")?no_esc}
    <#elseif section = "header">
        ${msg("loginProfileTitle")?no_esc}
    <#elseif section = "form">
        <form id="kc-update-profile-form" class="form update-profile ${properties.kcFormClass!}" action="${url.loginAction}" method="post">
            <#if user.editUsernameAllowed>
                <div class="update-profile-field ${properties.kcFormGroupClass!} ${messagesPerField.printIfExists('username',properties.kcFormGroupErrorClass!)}">
                    <div class="${properties.kcLabelWrapperClass!}">
                        <label for="username" class="${properties.kcLabelClass!}">${msg("username")}</label>
                    </div>
                    <div class="${properties.kcInputWrapperClass!}">
                        <input type="text" id="username" name="username" value="${(user.username!'')}" class="form-control ${properties.kcInputClass!}"/>
                    </div>
                </div>
            </#if>
            <div class="update-profile-field ${properties.kcFormGroupClass!} ${messagesPerField.printIfExists('email',properties.kcFormGroupErrorClass!)}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="email" class="${properties.kcLabelClass!}">${msg("email")}</label>
                </div>
                <div class="${properties.kcInputWrapperClass!}">
                    <input type="text" id="email" name="email" value="${(user.email!'')}" class="form-control ${properties.kcInputClass!}" />
                </div>
            </div>

            <div class="update-profile-field ${properties.kcFormGroupClass!} ${messagesPerField.printIfExists('firstName',properties.kcFormGroupErrorClass!)}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="firstName" class="${properties.kcLabelClass!}">${msg("firstName")}</label>
                </div>
                <div class="${properties.kcInputWrapperClass!}">
                    <input type="text" id="firstName" name="firstName" value="${(user.firstName!'')}" class="form-control ${properties.kcInputClass!}" />
                </div>
            </div>

            <div class="update-profile-field ${properties.kcFormGroupClass!} ${messagesPerField.printIfExists('lastName',properties.kcFormGroupErrorClass!)}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="lastName" class="${properties.kcLabelClass!}">${msg("lastName")}</label>
                </div>
                <div class="${properties.kcInputWrapperClass!}">
                    <input type="text" id="lastName" name="lastName" value="${(user.lastName!'')}" class="form-control ${properties.kcInputClass!}" />
                </div>
            </div>

            <div class="${properties.kcFormGroupClass!} row update-profile-button-container">
                <div id="kc-form-options" class="${properties.kcFormOptionsClass!}">
                    <div class="${properties.kcFormOptionsWrapperClass!}">
                    </div>
                </div>

                <div id="kc-form-buttons" class="${properties.kcFormButtonsClass!} col-xs-6 xs-xs-offset-6 col-sm-4 col-sm-offset-8">
                    <input class="btn btn-primary btn-flat btn-block ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" type="submit" value="${msg("doSubmit")}" />
                </div>
            </div>
        </form>
    </#if>
</@layout.registrationLayout>