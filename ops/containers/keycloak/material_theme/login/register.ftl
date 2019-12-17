<#import "template.ftl" as layout>
<@layout.registrationLayout; section>
    <#if section = "title">
        ${msg("Register a NOMAD account",(realm.displayName!''))?no_esc}
    <#elseif section = "header">
        <div class="title">
            ${msg("Register a NOMAD account",(realm.displayNameHtml!''))?no_esc}
        </div>
    <#elseif section = "form">
        <form id="kc-register-form" class="register form ${properties.kcFormClass!}" action="${url.registrationAction}" method="post">
          <input type="text" readonly value="this is not a login form" style="display: none;">
          <input type="password" readonly value="this is not a login form" style="display: none;">

          <#if !realm.registrationEmailAsUsername>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">person</i>
                <input id="username" class="mdc-text-field__input ${properties.kcInputClass!}" name="username" type="text" autofocus value="${(register.formData.username!'')}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="username" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("username")?no_esc}</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

          </#if>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">person</i>
                <input required id="firstName" class="mdc-text-field__input ${properties.kcInputClass!}" name="firstName" type="text" autofocus value="${(register.formData.firstName!'')}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="firstName" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("firstName")?no_esc}</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">person</i>
                <input required id="lastName" class="mdc-text-field__input ${properties.kcInputClass!}" name="lastName" type="text" value="${(register.formData.lastName!'')}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="lastName" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("lastName")?no_esc}</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">email</i>
                <input required id="email" class="mdc-text-field__input ${properties.kcInputClass!}" name="email" type="text" value="${(register.formData.email!'')}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="email" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("email")?no_esc}</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">room</i>
                <input required id="user.attributes.affiliation" class="mdc-text-field__input ${properties.kcInputClass!}" name="user.attributes.affiliation" type="text" value="${(register.formData.affiliation!'')}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="user.attributes.affiliation" class="mdc-floating-label ${properties.kcLabelClass!}">Affiliation</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">account_balance</i>
                <input required id="user.attributes.affiliation_address" class="mdc-text-field__input ${properties.kcInputClass!}" name="user.attributes.affiliation_address" type="text" value="${(register.formData.affiliation_address!'')}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="user.attributes.affiliation_address" class="mdc-floating-label ${properties.kcLabelClass!}">Affiliation address</label>
                </div>
                <div class="mdc-notched-outline">
                    <svg>
                        <path class="mdc-notched-outline__path"/>
                    </svg>
                </div>
                <div class="mdc-notched-outline__idle"></div>
            </div>

            <#if passwordRequired>

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">lock</i>
                <input required id="password" class="mdc-text-field__input ${properties.kcInputClass!}" name="password" type="password">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="password" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("password")?no_esc}</label>
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

            <#--  <div class="register-field ${properties.kcFormGroupClass!} ${messagesPerField.printIfExists('password-confirm',properties.kcFormGroupErrorClass!)}">
                <div class="${properties.kcLabelWrapperClass!}">
                    <label for="password-confirm" class="${properties.kcLabelClass!}">${msg("passwordConfirm")} *</label>
                </div>
                <div class="${properties.kcInputWrapperClass!}">
                    <input type="password" id="password-confirm" class="form-control ${properties.kcInputClass!}" name="password-confirm" />

                    <div id="capsLockWarning2" style="font-weight: bold; color: maroon; margin: 0 0 10px 0; display: none;">
                        <i class="fa fa-exclamation-triangle" style="color: #f0ad4e"></i>
                        ${msg("capsLockWarning")}
                    </div>
                </div>
            </div>  -->
            </#if>

            <#if recaptchaRequired??>
            <div class="form-group">
                <div class="${properties.kcInputWrapperClass!}">
                    <div class="g-recaptcha" data-size="compact" data-sitekey="${recaptchaSiteKey}"></div>
                </div>
            </div>
            </#if>

            <div class="${properties.kcFormGroupClass!} register-button-container">
                <div id="kc-form-buttons" class="${properties.kcFormButtonsClass!}">
                    <button class="mdc-button mdc-button--raised ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" type="submit">
                        ${msg("doRegister")?no_esc}
                    </button>
                </div>
            </div>
        </form>
    </#if>
</@layout.registrationLayout>