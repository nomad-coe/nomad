<#import "template.ftl" as layout>
<@layout.registrationLayout; section>

    <#if section = "title">
        ${msg("Edit your NOMAD account",(realm.displayName!''))?no_esc}
    <#elseif section = "header">
        <div class="title">
            ${msg("Edit your NOMAD account",(realm.displayNameHtml!''))?no_esc}
        </div>
    <#elseif section = "form">
        <form id="kc-register-form" class="register form ${properties.kcFormClass!}"  action="${url.accountUrl}" class="form-horizontal" method="post">
            <input type="text" readonly value="this is not a login form" style="display: none;">
            <input type="password" readonly value="this is not a login form" style="display: none;">

            <input type="hidden" id="stateChecker" name="stateChecker" value="${stateChecker}">

            <input type="hidden" id="username" name="username" value="${(account.username!'')}">

            <div class="mdc-text-field mdc-text-field--outlined mdc-text-field--with-leading-icon ${properties.kcLabelClass!}">
                <i class="material-icons mdc-text-field__icon" tabindex="-1" role="button">person</i>
                <input required id="firstName" class="mdc-text-field__input ${properties.kcInputClass!}" name="firstName" type="text" autofocus value="${(account.firstName!'')}">
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
                <input required id="lastName" class="mdc-text-field__input ${properties.kcInputClass!}" name="lastName" type="text" value="${(account.lastName!'')}">
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
                <input required id="email" class="mdc-text-field__input ${properties.kcInputClass!}" name="email" type="text" value="${(account.email!'')}">
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
                <input required id="user.attributes.affiliation" class="mdc-text-field__input ${properties.kcInputClass!}" name="user.attributes.affiliation" type="text" value="${(account.attributes.affiliation!'')}">
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
                <input required id="user.attributes.affiliation_address" class="mdc-text-field__input ${properties.kcInputClass!}" name="user.attributes.affiliation_address" type="text" value="${(account.attributes.affiliation_address!'')}">
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

            <div class="${properties.kcFormGroupClass!} register-button-container">
                <div id="kc-form-buttons" class="${properties.kcFormButtonsClass!}">
                    <button class="mdc-button mdc-button--raised ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" name="submitAction" type="submit" value="Save">
                        ${msg("doSave")?no_esc}
                    </button>
                </div>
            </div>
        </form>
    </#if>

</@layout.registrationLayout>