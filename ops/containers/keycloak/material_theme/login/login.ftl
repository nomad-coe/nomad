<#import "template.ftl" as layout>
<@layout.registrationLayout displayInfo=social.displayInfo; section>

    <#if section = "title">
        NOMAD login
    <#elseif section = "header">
	<div class="title">
        NOMAD login
    </div>
    <#elseif section = "form">
        <#if realm.password>
            <form id="kc-form-login" class="form ${properties.kcFormClass!}" action="${url.loginAction}" method="post">
                <div class="username-container ${properties.kcFormGroupClass!}">
                    <div class="${properties.kcInputWrapperClass!}">
                        <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!} <#if usernameEditDisabled??>mdc-text-field--disabled</#if>">
                            <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">person</i>
                            <input required id="username" class="mdc-text-field__input ${properties.kcInputClass!}" name="username" value="${(login.username!'')}" type="text" autofocus autocomplete="off" <#if usernameEditDisabled??>disabled</#if>>
                            <div class="${properties.kcLabelWrapperClass!}">
                                <label for="username" class="mdc-floating-label ${properties.kcLabelClass!}">
                                    <#if !realm.loginWithEmailAllowed>
                                        ${msg("username")}
                                    <#elseif !realm.registrationEmailAsUsername>
                                        ${msg("usernameOrEmail")}
                                    <#else>
                                        ${msg("email")}
                                    </#if>
                                </label>
                            </div>
                            <div class="mdc-notched-outline">
                                <svg>
                                    <path class="mdc-notched-outline__path"/>
                                </svg>
                            </div>
                            <div class="mdc-notched-outline__idle"></div>
                        </div>
                    </div>
                </div>
                <div class="password-container ${properties.kcFormGroupClass!}">
                    <div class="${properties.kcInputWrapperClass!}">

                        <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                            <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">lock</i>
                            <input required id="password" class="mdc-text-field__input ${properties.kcInputClass!}" name="password" type="password" autocomplete="off">
                            <div class="${properties.kcLabelWrapperClass!}">
                                <label for="password" class="mdc-floating-label ${properties.kcLabelClass!}">${msg("password")}</label>
                            </div>
                            <div class="mdc-notched-outline">
                                <svg>
                                    <path class="mdc-notched-outline__path"/>
                                </svg>
                            </div>
                            <div class="mdc-notched-outline__idle"></div>
                        </div>
                    </div>
                    <#--  <div id="capsLockWarning" style="font-weight: bold; color: maroon; margin: 0 0 10px 0; display: none;">
                      <i class="fa fa-exclamation-triangle" style="color: #f0ad4e"></i>
					  ${msg("capsLockWarning")}
                    </div>  -->
                </div>

                <div class="${properties.kcFormGroupClass!}">
                    <div id="kc-form-options" class="${properties.kcFormOptionsClass!}">
                        <div class="col-xs-7">
                            <#if realm.rememberMe && !usernameEditDisabled??>

                                <div class="mdc-form-field remember-me-checkbox">
                                    <div class="mdc-checkbox">
                                        <input type="checkbox"
                                               name="rememberMe"
                                               class="mdc-checkbox__native-control"
                                               id="rememberMe"
                                               <#if login.rememberMe??>checked</#if>
                                               />
                                        <div class="mdc-checkbox__background">
                                            <svg class="mdc-checkbox__checkmark"
                                                viewBox="0 0 24 24">
                                                <path class="mdc-checkbox__checkmark-path"
                                                    fill="none"
                                                    stroke="white"
                                                    d="M1.73,12.91 8.1,19.28 22.79,4.59"/>
                                            </svg>
                                            <div class="mdc-checkbox__mixedmark"></div>
                                        </div>
                                    </div>
                                    <label for="rememberMe">${msg("rememberMe")}</label>
                                </div>
                            </#if>
                        </div>
                        <div id="kc-form-buttons" class="col-xs-5 ${properties.kcFormButtonsClass!}">
                            <div class="${properties.kcFormButtonsWrapperClass!}">
                                <button class="mdc-button mdc-button--raised ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" name="login" id="kc-login" type="submit">
                                    ${msg("doLogIn")}
                                </button>
                            </div>
                            <div class="clearfix"></div>
                        </div>
                    </div>
                    <#if (realm.password && realm.registrationAllowed && !usernameEditDisabled??) || realm.resetPasswordAllowed>
                        <div>
                            <div class="col-xs-12">
                                <hr class="separator" />
                            </div>
                        </div>
                        <div>
                            <#if realm.password && realm.registrationAllowed && !usernameEditDisabled??>
                                <div id="kc-registration" class="col-xs-12" <#if realm.password && realm.registrationAllowed && !usernameEditDisabled?? && realm.resetPasswordAllowed>style="margin-bottom: 15px;"</#if>>
                                    <span>${msg("noAccount")} <a href="${url.registrationUrl}">${msg("doRegister")}</a></span>
                                </div>
                            </#if>
                        </div>
                        <div>
                            <div class="${properties.kcFormOptionsWrapperClass!} col-xs-12">
                                <#if realm.resetPasswordAllowed>
                                    <span><a href="${url.loginResetCredentialsUrl}">${msg("doForgotPassword")}</a></span>
                                </#if>
                            </div>
                        </div>
                    </#if>
                </div>
            </form>
        </#if>
    <#elseif section = "info" >
        <#if realm.password && social.providers??>
            <div id="kc-social-providers">
                <ul>
                    <#list social.providers as p>
                        <li><a href="${p.loginUrl}" id="zocial-${p.alias}" class="zocial ${p.providerId}"> <span class="text">${p.displayName}</span></a></li>
                    </#list>
                </ul>
            </div>
        </#if>
    </#if>
</@layout.registrationLayout>
